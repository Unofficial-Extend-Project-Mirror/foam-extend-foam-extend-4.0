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

Class
    Foam::Pstream

Description
    Inter-processor communications stream

SourceFiles
    Pstream.C
    PstreamCommsStruct.C
    gatherScatter.C
    combineGatherScatter.C
    gatherScatterList.C
    PstreamExchange.C

\*---------------------------------------------------------------------------*/

#ifndef Pstream_H
#define Pstream_H

#include "labelList.H"
#include "DynamicList.H"
#include "HashTable.H"
#include "foamString.H"
#include "NamedEnum.H"
#include "dynamicLabelList.H"
#include "optimisationSwitch.H"
#include "ListOps.H"
#include "LIFOStack.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Pstream Declaration
\*---------------------------------------------------------------------------*/

class Pstream
{
public:

    //- Types of communications
    enum commsTypes
    {
        blocking,
        scheduled,
        nonBlocking
    };

    static const NamedEnum<commsTypes, 3> commsTypeNames;


    // Public classes

        //- Structure for communicating between processors
        class commsStruct
        {
            // Private data

                //- procID of above processor
                label above_;

                //- procIDs of processors directly below me
                labelList below_;

                //- procIDs of all processors below (not just directly below)
                labelList allBelow_;

                //- procIDs of all processors not below. (inverse set of
                //  allBelow_ and minus myProcNo)
                labelList allNotBelow_;


        public:

            // Constructors

                //- Construct null
                commsStruct();

                //- Construct from components
                commsStruct
                (
                    const label,
                    const labelList&,
                    const labelList&,
                    const labelList&
                );

                //- Construct from components; construct allNotBelow_
                commsStruct
                (
                    const label nProcs,
                    const label myProcID,
                    const label,
                    const labelList&,
                    const labelList&
                );


            // Member Functions

                // Access

                    label above() const
                    {
                        return above_;
                    }

                    const labelList& below() const
                    {
                        return below_;
                    }

                    const labelList& allBelow() const
                    {
                        return allBelow_;
                    }

                    const labelList& allNotBelow() const
                    {
                        return allNotBelow_;
                    }


            // Member operators

                bool operator==(const commsStruct&) const;

                bool operator!=(const commsStruct&) const;


             // Ostream Operator

                friend Ostream& operator<<(Ostream&, const commsStruct&);
        };

        //- combineReduce operator for lists.  Used for counting.
        class listEq
        {
        public:

            template<class T>
            void operator()(T& x, const T& y) const
            {
                forAll (y, i)
                {
                    if (y[i].size())
                    {
                        x[i] = y[i];
                    }
                }
            }
        };


private:

    // Private data

        //- Is this a parallel run?
        static bool parRun_;

        //- Default message type info
        static const int msgType_;

        //- Stack of free comms
        static LIFOStack<label> freeComms_;

        //- My processor index
        static DynamicList<int> myProcNo_;

        //- Process IDs
        static DynamicList<List<int> > procIDs_;

        //- List of parent communicators
        static DynamicList<label> parentCommunicator_;

        //- Structure for linear communications
        static DynamicList<List<commsStruct> > linearCommunication_;

        //- Structure for tree communications
        static DynamicList<List<commsStruct> > treeCommunication_;


    // Private member functions

        //- Set data for parallel running
        static void setParRun(const label nProcs);

        //- Calculate linear communication schedule
        static List<commsStruct> calcLinearComm(const label nProcs);

        //- Calculate tree communication schedule
        static List<commsStruct> calcTreeComm(const label nProcs);

        //- Helper function for tree communication schedule determination
        //  Collects all processorIDs below a processor
        static void collectReceives
        (
            const label procID,
            const dynamicLabelListList& receives,
            dynamicLabelList& allReceives
        );

        //- Initialize all communication schedules. Callback from
        //  Pstream::init()
        static void initCommunicationSchedule();

        //- Allocate a communicator with index
        static void allocatePstreamCommunicator
        (
            const label parentIndex,
            const label index
        );

        //- Free a communicator
        static void freePstreamCommunicator
        (
            const label index
        );


protected:

    // Protected data

        //- Communications type of this stream
        commsTypes commsType_;

        //- Transfer buffer
        List<char> buf_;

        //- Current buffer read/write location
        int bufPosition_;


    // Protected member functions

        //- Increase the size of the transfer buffer
        inline void enlargeBuffer(size_t count);


public:

    // Declare name of the class and its debug switch
    ClassName("Pstream");


    // Static data

        //- Number of processors at which the sum algorithm changes from linear
        //  to tree
        static const debug::optimisationSwitch nProcsSimpleSum;

        //- Default commsType
        static debug::optimisationSwitch defaultCommsType;

        //- Default comms with a cast
        inline static Pstream::commsTypes defaultComms()
        {
            return static_cast<Pstream::commsTypes>
            (
                Pstream::defaultCommsType()
            );
        }

        //- Number of polling cycles in processor updates
        static const debug::optimisationSwitch nPollProcInterfaces;

        //- Default communicator (all processors)
        static label worldComm;

        //- Debugging: warn for use of any communicator differing from warnComm
        static label warnComm;


    // Constructors

        //- Construct given optional buffer size
        Pstream
        (
            const commsTypes commsType,
            const label bufSize = 0
        )
        :
            commsType_(commsType),
            bufPosition_(0)
        {
            if (bufSize)
            {
                buf_.setSize(bufSize + 2*sizeof(scalar) + 1);
            }
        }


    // Member functions

        //- Allocate a new communicator
        static label allocateCommunicator
        (
            const label parent,
            const labelList& subRanks,
            const bool doPstream = true
        );

        //- Free a previously allocated communicator
        static void freeCommunicator
        (
            const label communicator,
            const bool doPstream = true
        );

        //- Free all communicators
        static void freeCommunicators(const bool doPstream);


        //- Helper class for allocating/freeing communicators
        class communicator
        {
            //- Communicator identifier
            label comm_;

            //- Disallow copy and assignment
            communicator(const communicator&);
            void operator=(const communicator&);

        public:

            //- Constructo from components
            communicator
            (
                const label parent,
                const labelList& subRanks,
                const bool doPstream
            )
            :
                comm_(allocateCommunicator(parent, subRanks, doPstream))
            {}

            //- Destructor
            ~communicator()
            {
                freeCommunicator(comm_);
            }

            //- Cast to label
            operator label() const
            {
                return comm_;
            }
        };


        //- Return physical processor number (i.e. processor number in
        //  worldComm) given communicator and processor
        static int baseProcNo(const label myComm, const int procID);

        //- Return processor number in communicator (given physical processor
        //  number) (= reverse of baseProcNo)
        static label procNo(const label comm, const int baseProcID);

        //- Return processor number in communicator (given processor number
        //  and communicator)
        static label procNo
        (
            const label myComm,
            const label currentComm,
            const int currentProcID
        );

        //- Add the valid option this type of communications library
        //  adds/requires on the command line
        static void addValidParOptions(HashTable<string>& validParOptions);

        //- Initialisation function called from main
        //  Spawns slave processes and initialises inter-communication
        static bool init(int& argc, char**& argv);


        // Non-blocking comms

            //- Get number of outstanding requests
            static label nRequests();

            //- Truncate number of outstanding requests
            static void resetRequests(const label sz);

            //- Wait until all requests (from start onwards) have finished.
            static void waitRequests(const label start = 0);

            //- Wait until request i has finished.
            static void waitRequest(const label i);

            //- Non-blocking comms: has request i finished?
            static bool finishedRequest(const label i);

            static int allocateTag(const char*);

            static int allocateTag(const word&);

            static void freeTag(const char*, const int tag);

            static void freeTag(const word&, const int tag);


        //- Is this a parallel run?
        static bool& parRun()
        {
            return parRun_;
        }

        //- Number of processes in parallel run for a given communicator
        static label nProcs(const label communicator = 0)
        {
            return procIDs_[communicator].size();
        }

        //- Process index of the master for the global communicator
        static int masterNo()
        {
            return 0;
        }

        //- Am I the master process
        static bool master(const label communicator = 0)
        {
            return myProcNo_[communicator] == masterNo();
        }


        //- Number of this process (starting from masterNo() = 0)
        static int myProcNo(const label communicator = 0)
        {
            return myProcNo_[communicator];
        }

        //- Return parent communicator
        static label parent(const label communicator)
        {
            return parentCommunicator_(communicator);
        }

        //- Process IDs
        static const List<int>& procIDs(const int communicator)
        {
            return procIDs_[communicator];
        }

        //- Process ID of given process index
        static List<int>& procID(const int procNo)
        {
            return procIDs_[procNo];
        }

        //- Process index of first slave
        static int firstSlave()
        {
            return 1;
        }

        //- Process index of last slave
        static int lastSlave(const label communicator = 0)
        {
            return nProcs(communicator) - 1;
        }

        //- Communication schedule for linear all-to-master (proc 0)
        static const List<commsStruct>& linearCommunication
        (
            const label communicator = 0
        )
        {
            return linearCommunication_[communicator];
        }

        //- Communication schedule for tree all-to-master (proc 0)
        static const List<commsStruct>& treeCommunication
        (
            const label communicator = 0
        )
        {
            return treeCommunication_[communicator];
        }

        //- Message tag of standard messages
        static int msgType()
        {
            return msgType_;
        }

        //- Get the communications type of the stream
        commsTypes commsType() const
        {
            return commsType_;
        }

        //- Set the communications type of the stream
        commsTypes commsType(const commsTypes ct)
        {
            commsTypes oldCommsType = commsType_;
            commsType_ = ct;
            return oldCommsType;
        }

        //- Exit program
        static void exit(int errnum = 1);

        //- Abort program
        static void abort();


        // Gather and scatter

            //- Gather data. Apply bop to combine Value
            //  from different processors
            template <class T, class BinaryOp>
            static void gather
            (
                const List<commsStruct>& comms,
                T& Value,
                const BinaryOp& bop,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template <class T, class BinaryOp>
            static void gather
            (
                T& Value,
                const BinaryOp& bop,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );

            //- Scatter data. Distribute without modification.
            //  Reverse of gather
            template <class T>
            static void scatter
            (
                const List<commsStruct>& comms,
                T& Value,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template <class T>
            static void scatter
            (
                T& Value,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );


        // Combine variants. Inplace combine values from processors.
        // (Uses construct from Istream instead of <<)

            template <class T, class CombineOp>
            static void combineGather
            (
                const List<commsStruct>& comms,
                T& Value,
                const CombineOp& cop,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template <class T, class CombineOp>
            static void combineGather
            (
                T& Value,
                const CombineOp& cop,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );

            //- Scatter data. Reverse of combineGather
            template <class T>
            static void combineScatter
            (
                const List<commsStruct>& comms,
                T& Value,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template <class T>
            static void combineScatter
            (
                T& Value,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );

        // Combine variants working on whole List at a time.

            template <class T, class CombineOp>
            static void listCombineGather
            (
                const List<commsStruct>& comms,
                List<T>& Value,
                const CombineOp& cop,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template <class T, class CombineOp>
            static void listCombineGather
            (
                List<T>& Value,
                const CombineOp& cop,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );

            //- Scatter data. Reverse of combineGather
            template <class T>
            static void listCombineScatter
            (
                const List<commsStruct>& comms,
                List<T>& Value,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template <class T>
            static void listCombineScatter
            (
                List<T>& Value,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );

        // Combine variants working on whole map at a time. Container needs to
        // have iterators and find() defined.

            template <class Container, class CombineOp>
            static void mapCombineGather
            (
                const List<commsStruct>& comms,
                Container& Values,
                const CombineOp& cop,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template <class Container, class CombineOp>
            static void mapCombineGather
            (
                Container& Values,
                const CombineOp& cop,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );

            //- Scatter data. Reverse of combineGather
            template <class Container>
            static void mapCombineScatter
            (
                const List<commsStruct>& comms,
                Container& Values,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template <class Container>
            static void mapCombineScatter
            (
                Container& Values,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );


        // Gather/scatter keeping the individual processor data separate.
        // Values is a List of size Pstream::nProcs() where
        // Values[Pstream::myProcNo()] is the data for the current processor.

            //- Gather data but keep individual values separate
            template <class T>
            static void gatherList
            (
                const List<commsStruct>& comms,
                List<T>& Values,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template <class T>
            static void gatherList
            (
                List<T>& Values,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );

            //- Scatter data. Reverse of gatherList
            template <class T>
            static void scatterList
            (
                const List<commsStruct>& comms,
                List<T>& Values,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template <class T>
            static void scatterList
            (
                List<T>& Values,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );


        // Exchange

            //- Exchange data. Sends sendData, receives into recvData, sets
            //  sizes (not bytes). sizes[p0][p1] is what processor p0 has
            //  sent to p1. Continuous data only.
            //  If block=true will wait for all transfers to finish.
            template<class Container, class T>
            static void exchange
            (
                const List<Container>&,
                List<Container>&,
                labelListList& sizes,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm,
                const bool block = true
            );
};


inline void Pstream::enlargeBuffer(size_t count)
{
    buf_.setSize(max(int(buf_.size() + count), 2*buf_.size()));
}


Ostream& operator<<(Ostream&, const Pstream::commsStruct&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "gatherScatter.C"
#   include "combineGatherScatter.C"
#   include "gatherScatterList.C"
#   include "PstreamExchange.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
