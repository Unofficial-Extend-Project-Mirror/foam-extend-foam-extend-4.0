/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mpi.h"

#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "debug.H"
#include "dictionary.H"
#include "OSspecific.H"
#include "PstreamGlobals.H"
#include "SubList.H"

#include <cstring>
#include <cstdlib>
#include <csignal>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::Pstream, 0);

template<>
const char* Foam::NamedEnum<Foam::Pstream::commsTypes, 3>::names[] =
{
    "blocking",
    "scheduled",
    "nonBlocking"
};

const Foam::NamedEnum<Foam::Pstream::commsTypes, 3>
    Foam::Pstream::commsTypeNames;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Pstream::setParRun(const label nProcs)
{
    if (nProcs == 0)
    {
        parRun_ = false;
        freeCommunicator(Pstream::worldComm);

        label comm = allocateCommunicator(-1, labelList(1, label(0)), false);

        if (comm != Pstream::worldComm)
        {
            FatalErrorIn("Pstream::setParRun(const label)")
                << "problem : comm:" << comm
                << "  Pstream::worldComm:" << Pstream::worldComm
                << Foam::exit(FatalError);
        }

        Pout.prefix() = "";
        Perr.prefix() = "";
    }
    else
    {
        parRun_ = true;

        // Redo worldComm communicator (created at static initialisation)
        freeCommunicator(Pstream::worldComm);
        label comm = allocateCommunicator(-1, identity(nProcs), true);

        if (comm != Pstream::worldComm)
        {
            FatalErrorIn("Pstream::setParRun(const label)")
                << "problem : comm:" << comm
                << "  Pstream::worldComm:" << Pstream::worldComm
                << Foam::exit(FatalError);
        }

        Pout.prefix() = '[' +  name(myProcNo()) + "] ";
        Perr.prefix() = '[' +  name(myProcNo()) + "] ";
    }
}


Foam::List<Foam::Pstream::commsStruct> Foam::Pstream::calcLinearComm
(
    const label nProcs
)
{
    List<commsStruct> linearCommunication(nProcs);

    // Master
    labelList belowIDs(nProcs - 1);
    forAll (belowIDs, i)
    {
        belowIDs[i] = i + 1;
    }

    linearCommunication[0] = commsStruct
    (
        nProcs,
        0,
        -1,
        belowIDs,
        labelList()
    );

    // Slaves. Have no below processors, only communicate up to master
    for (label procID = 1; procID < nProcs; procID++)
    {
        linearCommunication[procID] = commsStruct
        (
            nProcs,
            procID,
            0,
            labelList(),
            labelList()
        );
    }

    return linearCommunication;
}


// Tree like schedule. For 8 procs:
// (level 0)
//      0 receives from 1
//      2 receives from 3
//      4 receives from 5
//      6 receives from 7
// (level 1)
//      0 receives from 2
//      4 receives from 6
// (level 2)
//      0 receives from 4
//
// The sends/receives for all levels are collected per processor (one send per
// processor; multiple receives possible) creating a table:
//
// So per processor:
// proc     receives from   sends to
// ----     -------------   --------
//  0       1,2,4           -
//  1       -               0
//  2       3               0
//  3       -               2
//  4       5               0
//  5       -               4
//  6       7               4
//  7       -               6
Foam::List<Foam::Pstream::commsStruct> Foam::Pstream::calcTreeComm
(
    const label nProcs
)
{
    label nLevels = 1;
    while ((1 << nLevels) < nProcs)
    {
        nLevels++;
    }

    List<dynamicLabelList> receives(nProcs);
    labelList sends(nProcs, -1);

    label offset = 2;
    label childOffset = offset/2;

    for (label level = 0; level < nLevels; level++)
    {
        label receiveID = 0;
        while (receiveID < nProcs)
        {
            // Determine processor that sends and we receive from
            label sendID = receiveID + childOffset;

            if (sendID < nProcs)
            {
                receives[receiveID].append(sendID);
                sends[sendID] = receiveID;
            }

            receiveID += offset;
        }

        offset <<= 1;
        childOffset <<= 1;
    }

    // For all processors find the processors it receives data from
    // (and the processors they receive data from etc.)
    List<dynamicLabelList> allReceives(nProcs);
    for (label procID = 0; procID < nProcs; procID++)
    {
        collectReceives(procID, receives, allReceives[procID]);
    }


    List<commsStruct> treeCommunication(nProcs);

    for (label procID = 0; procID < nProcs; procID++)
    {
        treeCommunication[procID] = commsStruct
        (
            nProcs,
            procID,
            sends[procID],
            receives[procID].shrink(),
            allReceives[procID].shrink()
        );
    }

    return treeCommunication;
}


// Append my children (and my children;s children etc.) to allReceives.
void Foam::Pstream::collectReceives
(
    const label procID,
    const List<dynamicLabelList>& receives,
    dynamicLabelList& allReceives
)
{
    const dynamicLabelList& myChildren = receives[procID];

    forAll (myChildren, childI)
    {
        allReceives.append(myChildren[childI]);
        collectReceives(myChildren[childI], receives, allReceives);
    }
}


// Callback from Pstream::init() : initialize linear and tree communication
// schedules now that nProcs is known.
void Foam::Pstream::initCommunicationSchedule()
{
    calcLinearComm(nProcs());
    calcTreeComm(nProcs());
}


void Foam::Pstream::allocatePstreamCommunicator
(
    const label parentIndex,
    const label index
)
{
    if (index == PstreamGlobals::MPIGroups_.size())
    {
        // Extend storage with dummy values
        MPI_Group newGroup;
        PstreamGlobals::MPIGroups_.append(newGroup);
        MPI_Comm newComm;
        PstreamGlobals::MPICommunicators_.append(newComm);
    }
    else if (index > PstreamGlobals::MPIGroups_.size())
    {
        FatalErrorIn
        (
            "Pstream::allocatePstreamCommunicator\n"
            "(\n"
            "    const label parentIndex,\n"
            "    const labelList& subRanks\n"
            ")\n"
        )   << "PstreamGlobals out of sync with Pstream data. Problem."
            << Foam::exit(FatalError);
    }


    if (parentIndex == -1)
    {
        // Allocate world communicator

        if (index != Pstream::worldComm)
        {
            FatalErrorIn
            (
                "Pstream::allocatePstreamCommunicator\n"
                "(\n"
                "    const label parentIndex,\n"
                "    const labelList& subRanks\n"
                ")\n"
            )   << "world communicator should always be index "
                << Pstream::worldComm << Foam::exit(FatalError);
        }

        PstreamGlobals::MPICommunicators_[index] = MPI_COMM_WORLD;
        MPI_Comm_group(MPI_COMM_WORLD, &PstreamGlobals::MPIGroups_[index]);
        MPI_Comm_rank
        (
            PstreamGlobals::MPICommunicators_[index],
            &myProcNo_[index]
        );

        // Set the number of processes to the actual number
        int numProcs;
        MPI_Comm_size(PstreamGlobals::MPICommunicators_[index], &numProcs);
        procIDs_[index] = identity(numProcs);
    }
    else
    {
        // Create new group
        MPI_Group_incl
        (
            PstreamGlobals::MPIGroups_[parentIndex],
            procIDs_[index].size(),
            procIDs_[index].begin(),
            &PstreamGlobals::MPIGroups_[index]
        );

        // Create new communicator
        MPI_Comm_create
        (
            PstreamGlobals::MPICommunicators_[parentIndex],
            PstreamGlobals::MPIGroups_[index],
            &PstreamGlobals::MPICommunicators_[index]
        );

        if (PstreamGlobals::MPICommunicators_[index] == MPI_COMM_NULL)
        {
            myProcNo_[index] = -1;
        }
        else
        {
            if
            (
                MPI_Comm_rank
                (
                    PstreamGlobals::MPICommunicators_[index],
                    &myProcNo_[index]
                )
            )
            {
                FatalErrorIn
                (
                    "Pstream::allocatePstreamCommunicator\n"
                    "(\n"
                    "    const label,\n"
                    "    const labelList&\n"
                    ")\n"
                )   << "Problem :"
                    << " when allocating communicator at " << index
                    << " from ranks " << procIDs_[index]
                    << " of parent " << parentIndex
                    << " cannot find my own rank"
                    << Foam::exit(FatalError);
            }
        }
    }
}


void Foam::Pstream::freePstreamCommunicator(const label communicator)
{
    if (communicator != Pstream::worldComm)
    {
        if (PstreamGlobals::MPICommunicators_[communicator] != MPI_COMM_NULL)
        {
            // Free communicator. Sets communicator to MPI_COMM_NULL
            MPI_Comm_free(&PstreamGlobals::MPICommunicators_[communicator]);
        }
        if (PstreamGlobals::MPIGroups_[communicator] != MPI_GROUP_NULL)
        {
            // Free greoup. Sets group to MPI_GROUP_NULL
            MPI_Group_free(&PstreamGlobals::MPIGroups_[communicator]);
        }
    }
}


Foam::label Foam::Pstream::allocateCommunicator
(
    const label parentIndex,
    const labelList& subRanks,
    const bool doPstream
)
{
    label index;
    if (!freeComms_.empty())
    {
        index = freeComms_.pop();
    }
    else
    {
        // Extend storage
        index = parentCommunicator_.size();

        myProcNo_.append(-1);
        procIDs_.append(List<int>());
        parentCommunicator_.append(-1);
        linearCommunication_.append(List<commsStruct>());
        treeCommunication_.append(List<commsStruct>());
    }

    if (debug)
    {
        Pout<< "Communicators : Allocating communicator " << index << endl
            << "    parent : " << parentIndex << endl
            << "    procs  : " << subRanks << endl
            << endl;
    }

    // Initialise; overwritten by allocatePstreamCommunicator
    myProcNo_[index] = 0;

    // Convert from label to int
    procIDs_[index].setSize(subRanks.size());
    forAll (procIDs_[index], i)
    {
        procIDs_[index][i] = subRanks[i];

        // Enforce incremental order (so index is rank in next communicator)
        if (i >= 1 && subRanks[i] <= subRanks[i - 1])
        {
            FatalErrorIn
            (
                "Pstream::allocateCommunicator"
                "(const label, const labelList&, const bool)"
            )   << "subranks not sorted : " << subRanks
                << " when allocating subcommunicator from parent "
                << parentIndex
                << Foam::abort(FatalError);
        }
    }
    parentCommunicator_[index] = parentIndex;

    linearCommunication_[index] = calcLinearComm(procIDs_[index].size());
    treeCommunication_[index] = calcTreeComm(procIDs_[index].size());


    if (doPstream && parRun())
    {
        allocatePstreamCommunicator(parentIndex, index);
    }

    return index;
}


void Foam::Pstream::freeCommunicator
(
    const label communicator,
    const bool doPstream
)
{
    if (debug)
    {
        Pout<< "Communicators : Freeing communicator " << communicator << endl
            << "    parent   : " << parentCommunicator_[communicator] << endl
            << "    myProcNo : " << myProcNo_[communicator] << endl
            << endl;
    }

    if (doPstream && parRun())
    {
        freePstreamCommunicator(communicator);
    }
    myProcNo_[communicator] = -1;
    //procIDs_[communicator].clear();
    parentCommunicator_[communicator] = -1;
    linearCommunication_[communicator].clear();
    treeCommunication_[communicator].clear();

    freeComms_.push(communicator);
}


void Foam::Pstream::freeCommunicators(const bool doPstream)
{
    forAll (myProcNo_, communicator)
    {
        if (myProcNo_[communicator] != -1)
        {
            freeCommunicator(communicator, doPstream);
        }
    }
}


int Foam::Pstream::baseProcNo(const label myComm, const int myProcID)
{
    int procID = myProcID;
    label comm = myComm;

    while (parent(comm) != -1)
    {
        const List<int>& parentRanks = Pstream::procID(comm);
        procID = parentRanks[procID];
        comm = Pstream::parent(comm);
    }

    return procID;
}


Foam::label Foam::Pstream::procNo(const label myComm, const int baseProcID)
{
    const List<int>& parentRanks = procID(myComm);
    label parentComm = parent(myComm);

    if (parentComm == -1)
    {
        return findIndex(parentRanks, baseProcID);
    }
    else
    {
        label parentRank = procNo(parentComm, baseProcID);
        return findIndex(parentRanks, parentRank);
    }
}


Foam::label Foam::Pstream::procNo
(
    const label myComm,
    const label currentComm,
    const int currentProcID
)
{
    label physProcID = Pstream::baseProcNo(currentComm, currentProcID);
    return procNo(myComm, physProcID);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// NOTE:
// valid parallel options vary between implementations, but flag common ones.
// if they are not removed by MPI_Init(), the subsequent argument processing
// will notice that they are wrong
void Foam::Pstream::addValidParOptions(HashTable<string>& validParOptions)
{
    validParOptions.insert("np", "");
    validParOptions.insert("p4pg", "PI file");
    validParOptions.insert("p4wd", "directory");
    validParOptions.insert("p4amslave", "");
    validParOptions.insert("p4yourname", "hostname");
    validParOptions.insert("machinefile", "machine file");
}


bool Foam::Pstream::init(int& argc, char**& argv)
{
    MPI_Init(&argc, &argv);

    int numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (debug)
    {
        Pout<< "Pstream::init : initialised with numProcs:" << numprocs
            << " myRank:" << myRank << endl;
    }

    if (numprocs <= 1)
    {
        FatalErrorIn("Pstream::init(int& argc, char**& argv)")
            << "bool IPstream::init(int& argc, char**& argv) : "
               "attempt to run parallel on 1 processor"
            << Foam::abort(FatalError);
    }

    // Initialise parallel structure
    setParRun(numprocs);

#   ifndef SGIMPI
    string bufferSizeName = getEnv("MPI_BUFFER_SIZE");

    if (bufferSizeName.size())
    {
        int bufferSize = atoi(bufferSizeName.c_str());

        if (bufferSize)
        {
            MPI_Buffer_attach(new char[bufferSize], bufferSize);
        }
    }
    else
    {
        FatalErrorIn("Pstream::init(int& argc, char**& argv)")
            << "Pstream::init(int& argc, char**& argv) : "
            << "environment variable MPI_BUFFER_SIZE not defined"
            << Foam::abort(FatalError);
    }
#   endif

    //int processorNameLen;
    //char processorName[MPI_MAX_PROCESSOR_NAME];
    //
    //MPI_Get_processor_name(processorName, &processorNameLen);
    //processorName[processorNameLen] = '\0';
    //Pout<< "Processor name:" << processorName << endl;

    return true;
}


void Foam::Pstream::exit(int errnum)
{
    if (debug)
    {
        Pout<< "Pstream::exit." << endl;
    }

#   ifndef SGIMPI
    int size;
    char* buff;
    MPI_Buffer_detach(&buff, &size);
    delete[] buff;
#   endif

    if (PstreamGlobals::outstandingRequests_.size())
    {
        label n = PstreamGlobals::outstandingRequests_.size();
        PstreamGlobals::outstandingRequests_.clear();

        WarningIn("Pstream::exit(int)")
            << "There are still " << n << " outstanding MPI_Requests." << endl
            << "This means that your code exited before doing a"
            << " Pstream::waitRequests()." << endl
            << "This should not happen for a normal code exit."
            << endl;
    }

    // Clean mpi communicators
    forAll (myProcNo_, communicator)
    {
        if (myProcNo_[communicator] != -1)
        {
            freePstreamCommunicator(communicator);
        }
    }

    if (errnum == 0)
    {
        MPI_Finalize();
        ::exit(errnum);
    }
    else
    {
        MPI_Abort(MPI_COMM_WORLD, errnum);
    }
}


void Foam::Pstream::abort()
{
    MPI_Abort(MPI_COMM_WORLD, 1);
}


Foam::label Foam::Pstream::nRequests()
{
    return PstreamGlobals::outstandingRequests_.size();
}


void Foam::Pstream::resetRequests(const label i)
{
    if (i < PstreamGlobals::outstandingRequests_.size())
    {
        PstreamGlobals::outstandingRequests_.setSize(i);
    }
}


void Foam::Pstream::waitRequests(const label start)
{
    if (debug)
    {
        Pout<< "Pstream::waitRequests : starting wait for "
            << PstreamGlobals::outstandingRequests_.size()-start
            << " outstanding requests starting at " << start << endl;
    }

    if (PstreamGlobals::outstandingRequests_.size())
    {
        SubList<MPI_Request> waitRequests
        (
            PstreamGlobals::outstandingRequests_,
            PstreamGlobals::outstandingRequests_.size() - start,
            start
        );

        if
        (
            MPI_Waitall
            (
                waitRequests.size(),
                waitRequests.begin(),
                MPI_STATUSES_IGNORE
            )
        )
        {
            FatalErrorIn
            (
                "Pstream::waitRequests()"
            )   << "MPI_Waitall returned with error" << Foam::endl;
        }

        resetRequests(start);
    }

    if (debug)
    {
        Pout<< "Pstream::waitRequests : finished wait." << endl;
    }
}


void Foam::Pstream::waitRequest(const label i)
{
    if (debug)
    {
        Pout<< "Pstream::waitRequest : starting wait for request:" << i
            << endl;
    }

    if (i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorIn
        (
            "Pstream::waitRequest(const label)"
        )   << "There are " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding send requests and you are asking for i=" << i
            << nl
            << "Maybe you are mixing blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }

    if
    (
        MPI_Wait
        (
            &PstreamGlobals::outstandingRequests_[i],
            MPI_STATUS_IGNORE
        )
    )
    {
        FatalErrorIn
        (
            "Pstream::waitRequest()"
        )   << "MPI_Wait returned with error" << Foam::endl;
    }

    if (debug)
    {
        Pout<< "Pstream::waitRequest : finished wait for request:" << i
            << endl;
    }
}


bool Foam::Pstream::finishedRequest(const label i)
{
    if (debug)
    {
        Pout<< "Pstream::finishedRequest : checking request:" << i
            << endl;
    }

    if (i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorIn
        (
            "Pstream::finishedRequest(const label)"
        )   << "There are " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding send requests and you are asking for i=" << i
            << nl
            << "Maybe you are mixing blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }

    int flag;
    MPI_Test
    (
        &PstreamGlobals::outstandingRequests_[i],
        &flag,
        MPI_STATUS_IGNORE
    );

    if (debug)
    {
        Pout<< "Pstream::finishedRequest : finished request:" << i
            << endl;
    }

    return flag != 0;
}


int Foam::Pstream::allocateTag(const char* s)
{
    int tag;
    if (PstreamGlobals::freedTags_.size())
    {
        tag = PstreamGlobals::freedTags_.remove();
    }
    else
    {
        tag = PstreamGlobals::nTags_++;
    }

    if (debug)
    {
        Pout<< "Pstream::allocateTag " << s
            << " : tag:" << tag
            << endl;
    }

    return tag;
}


int Foam::Pstream::allocateTag(const word& s)
{
    int tag;
    if (PstreamGlobals::freedTags_.size())
    {
        tag = PstreamGlobals::freedTags_.remove();
    }
    else
    {
        tag = PstreamGlobals::nTags_++;
    }

    if (debug)
    {
        Pout<< "Pstream::allocateTag " << s
            << " : tag:" << tag
            << endl;
    }

    return tag;
}


void Foam::Pstream::freeTag(const char* s, const int tag)
{
    if (debug)
    {
        Pout<< "Pstream::freeTag " << s << " tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.append(tag);
}


void Foam::Pstream::freeTag(const word& s, const int tag)
{
    if (debug)
    {
        Pout<< "Pstream::freeTag " << s << " tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.append(tag);
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// By default this is not a parallel run
bool Foam::Pstream::parRun_(false);

// Free communicators
Foam::LIFOStack<Foam::label> Foam::Pstream::freeComms_;

// My processor number
Foam::DynamicList<int> Foam::Pstream::myProcNo_(10);

// List of process IDs
Foam::DynamicList<Foam::List<int> > Foam::Pstream::procIDs_(10);

// Parent communicator
Foam::DynamicList<Foam::label> Foam::Pstream::parentCommunicator_(10);

// Standard transfer message type
const int Foam::Pstream::msgType_(1);

// Linear communication schedule
Foam::DynamicList<Foam::List<Foam::Pstream::commsStruct> >
Foam::Pstream::linearCommunication_(10);

// Multi level communication schedule
Foam::DynamicList<Foam::List<Foam::Pstream::commsStruct> >
Foam::Pstream::treeCommunication_(10);


// Allocate a serial communicator. This gets overwritten in parallel mode
// (by Pstream::setParRun())
Foam::Pstream::communicator serialComm(-1, Foam::labelList(1, 0), false);

// Number of processors at which the reduce algorithm changes from linear to
// tree
const Foam::debug::optimisationSwitch
Foam::Pstream::nProcsSimpleSum
(
    "nProcsSimpleSum",
    0
);


Foam::debug::optimisationSwitch
Foam::Pstream::defaultCommsType
(
    "commsType",
    "nonBlocking",
    "blocking, nonBlocking, scheduled"
);


const Foam::debug::optimisationSwitch
Foam::Pstream::nPollProcInterfaces
(
    "nPollProcInterfaces",
    0
);

// Default communicator
Foam::label Foam::Pstream::worldComm(0);


// Warn for use of any communicator
Foam::label Foam::Pstream::warnComm(-1);


// ************************************************************************* //
