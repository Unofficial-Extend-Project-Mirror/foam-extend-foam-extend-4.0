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

Description
    Gather data from all processors onto single processor according to some
    communication schedule (usually linear-to-master or tree-to-master).
    The gathered data will be a list with element procID the data from
     processor procID. Before calling every processor should insert
    its value into Values[Pstream::myProcNo()].

    Note: after gather every processor only knows its own data and that of the
    processors below it. Only the 'master' of the communication schedule holds
    a fully filled List. Use scatter to distribute the data.

\*---------------------------------------------------------------------------*/

#include "IPstream.H"
#include "OPstream.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T>
void Pstream::gatherList
(
    const List<Pstream::commsStruct>& comms,
    List<T>& Values,
    const int tag,
    const label comm
)
{
    // Return if not active in this comm
    // HJ, 12/Sep/2016
    if (Pstream::myProcNo(comm) == -1)
    {
        return;
    }

    if (Pstream::nProcs(comm) > 1)
    {
        if (Values.size() != Pstream::nProcs(comm))
        {
            FatalErrorIn
            (
                "void Pstream::gatherList\n"
                "(\n"
                "    const List<Pstream::commsStruct>& comms,\n"
                "    List<T>& Values,\n"
                "    const int tag,\n"
                "    const label comm\n"
                ")"
            )   << "Size of list:" << Values.size()
                << " does not equal the number of processors:"
                << Pstream::nProcs(comm)
                << Foam::abort(FatalError);
        }

        // Get my communication order
        const commsStruct& myComm = comms[Pstream::myProcNo(comm)];

        // Receive from my downstairs neighbours
        forAll (myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];
            const labelList& belowLeaves = comms[belowID].allBelow();

            if (contiguous<T>())
            {
                List<T> receivedValues(belowLeaves.size() + 1);

                IPstream::read
                (
                    Pstream::scheduled,
                    belowID,
                    reinterpret_cast<char*>(receivedValues.begin()),
                    receivedValues.byteSize(),
                    tag,
                    comm
                );

                Values[belowID] = receivedValues[0];

                forAll (belowLeaves, leafI)
                {
                    Values[belowLeaves[leafI]] = receivedValues[leafI + 1];
                }
            }
            else
            {
                IPstream fromBelow(Pstream::scheduled, belowID, 0, tag, comm);
                fromBelow >> Values[belowID];

                if (debug > 1)
                {
                    Pout<< " received through "
                        << belowID << " data from:" << belowID
                        << " data:" << Values[belowID] << endl;
                }

                // Receive from all other processors below belowID
                forAll (belowLeaves, leafI)
                {
                    label leafID = belowLeaves[leafI];
                    fromBelow >> Values[leafID];

                    if (debug > 1)
                    {
                        Pout<< " received through "
                            << belowID << " data from:" << leafID
                            << " data:" << Values[leafID] << endl;
                    }
                }
            }
        }

        // Send up from Values:
        // - my own value first
        // - all belowLeaves next
        if (myComm.above() != -1)
        {
            const labelList& belowLeaves = myComm.allBelow();

            if (debug > 1)
            {
                Pout<< " sending to " << myComm.above()
                    << " data from: " << Pstream::myProcNo(comm)
                    << " data: " << Values[Pstream::myProcNo(comm)] << endl;
            }

            if (contiguous<T>())
            {
                List<T> sendingValues(belowLeaves.size() + 1);
                sendingValues[0] = Values[Pstream::myProcNo(comm)];

                forAll (belowLeaves, leafI)
                {
                    sendingValues[leafI + 1] = Values[belowLeaves[leafI]];
                }

                OPstream::write
                (
                    Pstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<const char*>(sendingValues.begin()),
                    sendingValues.byteSize(),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream toAbove
                (
                    Pstream::scheduled,
                    myComm.above(),
                    0,
                    tag,
                    comm
                );
                toAbove << Values[Pstream::myProcNo(comm)];

                forAll (belowLeaves, leafI)
                {
                    label leafID = belowLeaves[leafI];

                    if (debug > 1)
                    {
                        Pout<< " sending to "
                            << myComm.above() << " data from: " << leafID
                            << " data: " << Values[leafID] << endl;
                    }
                    toAbove << Values[leafID];
                }
            }
        }
    }
}


template <class T>
void Pstream::gatherList(List<T>& Values, const int tag, const label comm)
{
    if (Pstream::nProcs(comm) < Pstream::nProcsSimpleSum())
    {
        gatherList(Pstream::linearCommunication(comm), Values, tag, comm);
    }
    else
    {
        gatherList(Pstream::treeCommunication(comm), Values, tag, comm);
    }
}


template <class T>
void Pstream::scatterList
(
    const List<Pstream::commsStruct>& comms,
    List<T>& Values,
    const int tag,
    const label comm
)
{
    // Return if not active in this comm
    // HJ, 12/Sep/2016
    if (Pstream::myProcNo(comm) == -1)
    {
        return;
    }

    if (Pstream::nProcs(comm) > 1)
    {
        if (Values.size() != Pstream::nProcs(comm))
        {
            FatalErrorIn
            (
                "void Pstream::scatterList\n"
                "(\n"
                "    const List<Pstream::commsStruct>& comms,\n"
                "    List<T>& Values,\n"
                "    const int tag,\n"
                "    const label comm\n"
                ")"
            )   << "Size of list:" << Values.size()
                << " does not equal the number of processors:"
                << Pstream::nProcs(comm)
                << Foam::abort(FatalError);
        }

        // Get my communication order
        const commsStruct& myComm = comms[Pstream::myProcNo(comm)];

        // Receive from up
        if (myComm.above() != -1)
        {
            const labelList& notBelowLeaves = myComm.allNotBelow();

            if (contiguous<T>())
            {
                List<T> receivedValues(notBelowLeaves.size());

                IPstream::read
                (
                    Pstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<char*>(receivedValues.begin()),
                    receivedValues.byteSize(),
                    tag,
                    comm
                );

                forAll (notBelowLeaves, leafI)
                {
                    Values[notBelowLeaves[leafI]] = receivedValues[leafI];
                }
            }
            else
            {
                IPstream fromAbove
                (
                    Pstream::scheduled,
                    myComm.above(),
                    0,
                    tag,
                    comm
                );

                forAll (notBelowLeaves, leafI)
                {
                    label leafID = notBelowLeaves[leafI];
                    fromAbove >> Values[leafID];

                    if (debug > 1)
                    {
                        Pout<< " received through "
                            << myComm.above() << " data for:" << leafID
                            << " data:" << Values[leafID] << endl;
                    }
                }
            }
        }

        // Send to my downstairs neighbours.  Note reverse order (compared to
        // receiving). This is to make sure to send to the critical path
        // (only when using a tree schedule!) first.
        // This is ESI Comms optimisation, v16.06.  HJ, 19/Sep/2016
        forAllReverse (myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];
            const labelList& notBelowLeaves = comms[belowID].allNotBelow();

            if (contiguous<T>())
            {
                List<T> sendingValues(notBelowLeaves.size());

                forAll (notBelowLeaves, leafI)
                {
                    sendingValues[leafI] = Values[notBelowLeaves[leafI]];
                }

                OPstream::write
                (
                    Pstream::scheduled,
                    belowID,
                    reinterpret_cast<const char*>(sendingValues.begin()),
                    sendingValues.byteSize(),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream toBelow(Pstream::scheduled, belowID, 0, tag, comm);

                // Send data destined for all other processors below belowID
                forAll (notBelowLeaves, leafI)
                {
                    label leafID = notBelowLeaves[leafI];
                    toBelow << Values[leafID];

                    if (debug > 1)
                    {
                        Pout<< " sent through "
                            << belowID << " data for:" << leafID
                            << " data:" << Values[leafID] << endl;
                    }
                }
            }
        }
    }
}


template <class T>
void Pstream::scatterList(List<T>& Values, const int tag, const label comm)
{
    if (Pstream::nProcs(comm) < Pstream::nProcsSimpleSum())
    {
        scatterList(Pstream::linearCommunication(comm), Values, tag, comm);
    }
    else
    {
        scatterList(Pstream::treeCommunication(comm), Values, tag, comm);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
