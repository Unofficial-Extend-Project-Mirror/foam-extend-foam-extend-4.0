/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "mapDistribute.H"
#include "commSchedule.H"
#include "HashSet.H"
#include "globalIndex.H"
#include "ListOps.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::mapDistribute, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::List<Foam::labelPair> Foam::mapDistribute::schedule
(
    const labelListList& subMap,
    const labelListList& constructMap,
    const int tag
)
{
    // Communications: send and receive processor
    List<labelPair> allComms;

    {
        HashSet<labelPair, labelPair::Hash<> > commsSet(Pstream::nProcs());

        // Find what communication is required
        forAll (subMap, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                if (subMap[procI].size())
                {
                    // I need to send to procI
                    commsSet.insert(labelPair(Pstream::myProcNo(), procI));
                }
                if (constructMap[procI].size())
                {
                    // I need to receive from procI
                    commsSet.insert(labelPair(procI, Pstream::myProcNo()));
                }
            }
        }
        allComms = commsSet.toc();
    }


    // Reduce
    if (Pstream::master())
    {
        // Receive and merge
        for
        (
            int slave = Pstream::firstSlave();
            slave <= Pstream::lastSlave();
            slave++
        )
        {
            IPstream fromSlave
            (
                Pstream::scheduled,
                slave,
                0,
                tag
            );

            List<labelPair> nbrData(fromSlave);

            forAll (nbrData, i)
            {
                if (findIndex(allComms, nbrData[i]) == -1)
                {
                    label sz = allComms.size();
                    allComms.setSize(sz + 1);
                    allComms[sz] = nbrData[i];
                }
            }
        }
        // Send back
        for
        (
            int slave = Pstream::firstSlave();
            slave <= Pstream::lastSlave();
            slave++
        )
        {
            OPstream toSlave
            (
                Pstream::scheduled,
                slave,
                0,
                tag
            );

            toSlave << allComms;
        }
    }
    else
    {
        {
            OPstream toMaster
            (
                Pstream::scheduled,
                Pstream::masterNo(),
                0,
                tag
            );

            toMaster << allComms;
        }
        {
            IPstream fromMaster
            (
                Pstream::scheduled,
                Pstream::masterNo(),
                0,
                tag
            );
            fromMaster >> allComms;
        }
    }


    // Determine my schedule.
    labelList mySchedule
    (
        commSchedule
        (
            Pstream::nProcs(),
            allComms
        ).procSchedule()[Pstream::myProcNo()]
    );

    // Processors involved in my schedule
    return List<labelPair>(UIndirectList<labelPair>(allComms, mySchedule));
}


const Foam::List<Foam::labelPair>& Foam::mapDistribute::schedule() const
{
    if (schedulePtr_.empty())
    {
        schedulePtr_.reset
        (
            new List<labelPair>
            (
                schedule(subMap_, constructMap_, Pstream::msgType())
            )
        );
    }

    return schedulePtr_();
}


void Foam::mapDistribute::checkReceivedSize
(
    const label procI,
    const label expectedSize,
    const label receivedSize
)
{
    if (receivedSize != expectedSize)
    {
        FatalErrorIn
        (
            "template<class T>\n"
            "void mapDistribute::checkReceivedSize\n"
            "(\n"
            "    const label procI,\n"
            "    const label expectedSize,\n"
            "    const label receivedSize\n"
            ")\n"
        )   << "Expected from processor " << procI
            << " " << expectedSize << " but received "
            << receivedSize << " elements."
            << abort(FatalError);
    }
}


// Construct per processor compact addressing of the global elements
// needed. The ones from the local processor are not included since
// these are always all needed.
void Foam::mapDistribute::calcCompactAddressing
(
    const globalIndex& globalNumbering,
    const labelList& elements,
    List<Map<label> >& compactMap
) const
{
    compactMap.setSize(Pstream::nProcs());

    // Count all (non-local) elements needed. Just for presizing map.
    labelList nNonLocal(Pstream::nProcs(), 0);

    forAll (elements, i)
    {
        label globalIndex = elements[i];

        if (globalIndex != -1 && !globalNumbering.isLocal(globalIndex))
        {
            label procI = globalNumbering.whichProcID(globalIndex);
            nNonLocal[procI]++;
        }
    }

    forAll (compactMap, procI)
    {
        compactMap[procI].clear();
        if (procI != Pstream::myProcNo())
        {
            compactMap[procI].resize(2*nNonLocal[procI]);
        }
    }


    // Collect all (non-local) elements needed.
    forAll (elements, i)
    {
        label globalIndex = elements[i];

        if (globalIndex != -1 && !globalNumbering.isLocal(globalIndex))
        {
            label procI = globalNumbering.whichProcID(globalIndex);
            label index = globalNumbering.toLocal(procI, globalIndex);
            label nCompact = compactMap[procI].size();
            compactMap[procI].insert(index, nCompact);
        }
    }
}


void Foam::mapDistribute::calcCompactAddressing
(
    const globalIndex& globalNumbering,
    const labelListList& cellCells,
    List<Map<label> >& compactMap
) const
{
    compactMap.setSize(Pstream::nProcs());

    // Count all (non-local) elements needed. Just for presizing map.
    labelList nNonLocal(Pstream::nProcs(), 0);

    forAll (cellCells, cellI)
    {
        const labelList& cCells = cellCells[cellI];

        forAll (cCells, i)
        {
            label globalIndex = cCells[i];

            if (globalIndex != -1 && !globalNumbering.isLocal(globalIndex))
            {
                label procI = globalNumbering.whichProcID(globalIndex);
                nNonLocal[procI]++;
            }
        }
    }

    forAll (compactMap, procI)
    {
        compactMap[procI].clear();
        if (procI != Pstream::myProcNo())
        {
            compactMap[procI].resize(2*nNonLocal[procI]);
        }
    }


    // Collect all (non-local) elements needed.
    forAll (cellCells, cellI)
    {
        const labelList& cCells = cellCells[cellI];

        forAll (cCells, i)
        {
            label globalIndex = cCells[i];

            if (globalIndex != -1 && !globalNumbering.isLocal(globalIndex))
            {
                label procI = globalNumbering.whichProcID(globalIndex);
                label index = globalNumbering.toLocal(procI, globalIndex);
                label nCompact = compactMap[procI].size();
                compactMap[procI].insert(index, nCompact);
            }
        }
    }
}


void Foam::mapDistribute::exchangeAddressing
(
    const int tag,
    const globalIndex& globalNumbering,
    labelList& elements,
    List<Map<label> >& compactMap,
    labelList& compactStart
)
{
    // The overall compact addressing is
    // - myProcNo data first (uncompacted)
    // - all other processors consecutively

    compactStart.setSize(Pstream::nProcs());
    compactStart[Pstream::myProcNo()] = 0;
    constructSize_ = globalNumbering.localSize();
    forAll (compactStart, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            compactStart[procI] = constructSize_;
            constructSize_ += compactMap[procI].size();
        }
    }

    // Find out what to receive/send in compact addressing.

    // What I want to receive is what others have to send
    labelListList wantedRemoteElements(Pstream::nProcs());
    // Compact addressing for received data
    constructMap_.setSize(Pstream::nProcs());
    forAll (compactMap, procI)
    {
        if (procI == Pstream::myProcNo())
        {
            // All my own elements are used
            label nLocal = globalNumbering.localSize();
            wantedRemoteElements[procI] = identity(nLocal);
            constructMap_[procI] = identity(nLocal);
        }
        else
        {
            // Remote elements wanted from processor procI
            labelList& remoteElem = wantedRemoteElements[procI];
            labelList& localElem = constructMap_[procI];
            remoteElem.setSize(compactMap[procI].size());
            localElem.setSize(compactMap[procI].size());
            label i = 0;
            forAllIter(Map<label>, compactMap[procI], iter)
            {
                const label compactI = compactStart[procI] + iter();
                remoteElem[i] = iter.key();
                localElem[i]  = compactI;
                iter() = compactI;
                i++;
            }
        }
    }

    subMap_.setSize(Pstream::nProcs());
    labelListList sendSizes;
    Pstream::exchange<labelList, label>
    (
        wantedRemoteElements,
        subMap_,
        sendSizes,
        tag
    );

    // Renumber elements
    forAll (elements, i)
    {
        elements[i] = renumber(globalNumbering, compactMap, elements[i]);
    }
}


void Foam::mapDistribute::exchangeAddressing
(
    const int tag,
    const globalIndex& globalNumbering,
    labelListList& cellCells,
    List<Map<label> >& compactMap,
    labelList& compactStart
)
{
    // The overall compact addressing is
    // - myProcNo data first (uncompacted)
    // - all other processors consecutively

    compactStart.setSize(Pstream::nProcs());
    compactStart[Pstream::myProcNo()] = 0;
    constructSize_ = globalNumbering.localSize();
    forAll (compactStart, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            compactStart[procI] = constructSize_;
            constructSize_ += compactMap[procI].size();
        }
    }

    // Find out what to receive/send in compact addressing.

    // What I want to receive is what others have to send
    labelListList wantedRemoteElements(Pstream::nProcs());
    // Compact addressing for received data
    constructMap_.setSize(Pstream::nProcs());
    forAll (compactMap, procI)
    {
        if (procI == Pstream::myProcNo())
        {
            // All my own elements are used
            label nLocal = globalNumbering.localSize();
            wantedRemoteElements[procI] = identity(nLocal);
            constructMap_[procI] = identity(nLocal);
        }
        else
        {
            // Remote elements wanted from processor procI
            labelList& remoteElem = wantedRemoteElements[procI];
            labelList& localElem = constructMap_[procI];
            remoteElem.setSize(compactMap[procI].size());
            localElem.setSize(compactMap[procI].size());
            label i = 0;
            forAllIter(Map<label>, compactMap[procI], iter)
            {
                const label compactI = compactStart[procI] + iter();
                remoteElem[i] = iter.key();
                localElem[i]  = compactI;
                iter() = compactI;
                i++;
            }
        }
    }

    subMap_.setSize(Pstream::nProcs());
    labelListList sendSizes;
    Pstream::exchange<labelList, label>
    (
        wantedRemoteElements,
        subMap_,
        sendSizes,
        tag
    );

    // Renumber elements
    forAll (cellCells, cellI)
    {
        labelList& cCells = cellCells[cellI];

        forAll (cCells, i)
        {
            cCells[i] = renumber(globalNumbering, compactMap, cCells[i]);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
Foam::mapDistribute::mapDistribute
(
    const label constructSize,
    const labelListList& subMap,
    const labelListList& constructMap
)
:
    constructSize_(constructSize),
    subMap_(subMap),
    constructMap_(constructMap),
    schedulePtr_()
{}


//- (optionally destructively) construct from components
Foam::mapDistribute::mapDistribute
(
    const label constructSize,
    labelListList& subMap,
    labelListList& constructMap,
    const bool reUse                // clone or reuse
)
:
    constructSize_(constructSize),
    subMap_(subMap, reUse),
    constructMap_(constructMap, reUse),
    schedulePtr_()
{}


Foam::mapDistribute::mapDistribute
(
    const labelList& sendProcs,
    const labelList& recvProcs
)
:
    constructSize_(0),
    schedulePtr_()
{
    if (sendProcs.size() != recvProcs.size())
    {
        FatalErrorIn
        (
            "mapDistribute::mapDistribute(const labelList&, const labelList&)"
        )   << "The send and receive data is not the same length. sendProcs:"
            << sendProcs.size() << " recvProcs:" << recvProcs.size()
            << abort(FatalError);
    }

    // Per processor the number of samples we have to send/receive.
    labelList nSend(Pstream::nProcs(), 0);
    labelList nRecv(Pstream::nProcs(), 0);

    forAll (sendProcs, sampleI)
    {
        label sendProc = sendProcs[sampleI];
        label recvProc = recvProcs[sampleI];

        // Note that also need to include local communication (both
        // RecvProc and sendProc on local processor)

        if (Pstream::myProcNo() == sendProc)
        {
            // I am the sender. Count destination processor.
            nSend[recvProc]++;
        }
        if (Pstream::myProcNo() == recvProc)
        {
            // I am the receiver.
            nRecv[sendProc]++;
        }
    }

    subMap_.setSize(Pstream::nProcs());
    constructMap_.setSize(Pstream::nProcs());

    forAll (nSend, procI)
    {
        subMap_[procI].setSize(nSend[procI]);
        constructMap_[procI].setSize(nRecv[procI]);
    }
    nSend = 0;
    nRecv = 0;

    forAll (sendProcs, sampleI)
    {
        label sendProc = sendProcs[sampleI];
        label recvProc = recvProcs[sampleI];

        if (Pstream::myProcNo() == sendProc)
        {
            // I am the sender. Store index I need to send.
            subMap_[recvProc][nSend[recvProc]++] = sampleI;
        }

        if (Pstream::myProcNo() == recvProc)
        {
            // I am the receiver.
            constructMap_[sendProc][nRecv[sendProc]++] = sampleI;

            // Largest entry inside constructMap
            constructSize_ = sampleI + 1;
        }
    }
}


Foam::mapDistribute::mapDistribute
(
    const globalIndex& globalNumbering,
    labelList& elements,
    List<Map<label> >& compactMap,
    const int tag
)
:
    constructSize_(0),
    schedulePtr_()
{
    // Construct per processor compact addressing of the global elements
    // needed. The ones from the local processor are not included since
    // these are always all needed.
    calcCompactAddressing
    (
        globalNumbering,
        elements,
        compactMap
    );

    // Exchange what I need with processor that supplies it. Renumber elements
    // into compact numbering
    labelList compactStart;
    exchangeAddressing
    (
        tag,
        globalNumbering,
        elements,
        compactMap,
        compactStart
    );

    if (debug)
    {
        printLayout(Pout);
    }
}


Foam::mapDistribute::mapDistribute
(
    const globalIndex& globalNumbering,
    labelListList& cellCells,
    List<Map<label> >& compactMap,
    const int tag
)
:
    constructSize_(0),
    schedulePtr_()
{
    // Construct per processor compact addressing of the global elements
    // needed. The ones from the local processor are not included since
    // these are always all needed.
    calcCompactAddressing
    (
        globalNumbering,
        cellCells,
        compactMap
    );

    // Exchange what I need with processor that supplies it. Renumber elements
    // into compact numbering
    labelList compactStart;

    exchangeAddressing
    (
        tag,
        globalNumbering,
        cellCells,
        compactMap,
        compactStart
    );

    if (debug)
    {
        printLayout(Pout);
    }
}


Foam::mapDistribute::mapDistribute(const mapDistribute& map)
:
    constructSize_(map.constructSize_),
    subMap_(map.subMap_),
    constructMap_(map.constructMap_),
    schedulePtr_()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::mapDistribute::renumber
(
    const globalIndex& globalNumbering,
    const List<Map<label> >& compactMap,
    const label globalI
)
{
    if (globalI == -1)
    {
        return globalI;
    }
    if (globalNumbering.isLocal(globalI))
    {
        return globalNumbering.toLocal(globalI);
    }
    else
    {
        label procI = globalNumbering.whichProcID(globalI);
        label index = globalNumbering.toLocal(procI, globalI);
        return compactMap[procI][index];
    }
}


void Foam::mapDistribute::compact(const boolList& elemIsUsed, const int tag)
{
    // 1. send back to sender. Have him delete the corresponding element
    //    from the submap and do the same to the constructMap locally
    //    (and in same order).

    // Send elemIsUsed field to neighbour. Use nonblocking code from
    // mapDistribute but in reverse order.
    if (Pstream::parRun())
    {
        label startOfRequests = Pstream::nRequests();

        // Set up receives from neighbours

        List<boolList> sendFields(Pstream::nProcs());

        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = constructMap_[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                boolList& subField = sendFields[domain];
                subField.setSize(map.size());
                forAll (map, i)
                {
                    subField[i] = elemIsUsed[map[i]];
                }

                OPstream::write
                (
                    Pstream::nonBlocking,
                    domain,
                    reinterpret_cast<const char*>(subField.begin()),
                    subField.size()*sizeof(bool),
                    tag
                );
            }
        }

        // Set up receives from neighbours

        List<boolList> recvFields(Pstream::nProcs());

        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = subMap_[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                recvFields[domain].setSize(map.size());
                IPstream::read
                (
                    Pstream::nonBlocking,
                    domain,
                    reinterpret_cast<char*>(recvFields[domain].begin()),
                    recvFields[domain].size()*sizeof(bool),
                    tag
                );
            }
        }


        // Set up 'send' to myself - write directly into recvFields

        {
            const labelList& map = constructMap_[Pstream::myProcNo()];

            recvFields[Pstream::myProcNo()].setSize(map.size());
            forAll (map, i)
            {
                recvFields[Pstream::myProcNo()][i] = elemIsUsed[map[i]];
            }
        }

        // Wait for all to finish
        Pstream::waitRequests(startOfRequests);

        // Compact out all submap entries that are referring to unused elements
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = subMap_[domain];

            labelList newMap(map.size());
            label newI = 0;

            forAll (map, i)
            {
                if (recvFields[domain][i])
                {
                    // So element is used on destination side
                    newMap[newI++] = map[i];
                }
            }
            if (newI < map.size())
            {
                newMap.setSize(newI);
                subMap_[domain].transfer(newMap);
            }
        }
    }


    // 2. remove from construct map - since end-result (element in elemIsUsed)
    //    not used.

    label maxConstructIndex = -1;

    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& map = constructMap_[domain];

        labelList newMap(map.size());
        label newI = 0;

        forAll (map, i)
        {
            label destinationI = map[i];

            // Is element is used on destination side
            if (elemIsUsed[destinationI])
            {
                maxConstructIndex = max(maxConstructIndex, destinationI);

                newMap[newI++] = destinationI;
            }
        }
        if (newI < map.size())
        {
            newMap.setSize(newI);
            constructMap_[domain].transfer(newMap);
        }
    }

    constructSize_ = maxConstructIndex + 1;

    // Clear the schedule (note:not necessary if nothing changed)
    Pout<< "CLEAR SCHEDULE POINTER" << endl;
    schedulePtr_.clear();
}


void Foam::mapDistribute::printLayout(Ostream& os) const
{
    // Determine offsets of remote data.
    labelList minIndex(Pstream::nProcs(), labelMax);
    labelList maxIndex(Pstream::nProcs(), labelMin);
    forAll (constructMap_, procI)
    {
        const labelList& construct = constructMap_[procI];
        minIndex[procI] = min(minIndex[procI], min(construct));
        maxIndex[procI] = max(maxIndex[procI], max(construct));
    }

    label localSize;
    if (maxIndex[Pstream::myProcNo()] == labelMin)
    {
        localSize = 0;
    }
    else
    {
        localSize = maxIndex[Pstream::myProcNo()]+1;
    }

    os  << "Layout: (constructSize:" << constructSize_ << ")" << endl
        << "local (processor " << Pstream::myProcNo() << "):" << endl
        << "    start : 0" << endl
        << "    size  : " << localSize << endl;

    label offset = localSize;
    forAll (minIndex, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            if (constructMap_[procI].size() > 0)
            {
                if (minIndex[procI] != offset)
                {
                    FatalErrorIn("mapDistribute::printLayout(..)")
                        << "offset:" << offset
                        << " procI:" << procI
                        << " minIndex:" << minIndex[procI]
                        << abort(FatalError);
                }

                label size = maxIndex[procI]-minIndex[procI]+1;
                os  << "processor " << procI << ':' << endl
                    << "    start : " << offset << endl
                    << "    size  : " << size << endl;

                offset += size;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::mapDistribute::operator=(const mapDistribute& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "Foam::mapDistribute::operator=(const Foam::mapDistribute&)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }
    constructSize_ = rhs.constructSize_;
    subMap_ = rhs.subMap_;
    constructMap_ = rhs.constructMap_;
    schedulePtr_.clear();
}


// ************************************************************************* //
