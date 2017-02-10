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
    Variant of gather, scatter.
    Normal gather uses:
    - construct null and read (>>) from Istream
    - binary operator and assignment operator to combine values

    combineGather uses:
    - construct from Istream
    - modify operator which modifies its lhs

\*---------------------------------------------------------------------------*/

#include "OPstream.H"
#include "IPstream.H"
#include "IOstreams.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T, class CombineOp>
void Pstream::combineGather
(
    const List<Pstream::commsStruct>& comms,
    T& Value,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (Pstream::nProcs(comm) > 1)
    {
        // Get my communication order
        const commsStruct& myComm = comms[Pstream::myProcNo(comm)];

        // Receive from my downstairs neighbours
        forAll (myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            if (contiguous<T>())
            {
                T value;
                IPstream::read
                (
                    Pstream::scheduled,
                    belowID,
                    reinterpret_cast<char*>(&value),
                    sizeof(T),
                    tag,
                    comm
                );

                if (debug > 1)
                {
                    Pout<< " received from "
                        << belowID << " data:" << value << endl;
                }

                cop(Value, value);
            }
            else
            {
                IPstream fromBelow(Pstream::scheduled, belowID, 0, tag, comm);
                T value(fromBelow);

                if (debug > 1)
                {
                    Pout<< " received from "
                        << belowID << " data:" << value << endl;
                }

                cop(Value, value);
            }
        }

        // Send up Value
        if (myComm.above() != -1)
        {
            if (debug > 1)
            {
                Pout<< " sending to " << myComm.above()
                    << " data:" << Value << endl;
            }

            if (contiguous<T>())
            {
                OPstream::write
                (
                    Pstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<const char*>(&Value),
                    sizeof(T),
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
                toAbove << Value;
            }
        }
    }
}


template <class T, class CombineOp>
void Pstream::combineGather
(
    T& Value,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (Pstream::nProcs(comm) < Pstream::nProcsSimpleSum())
    {
        combineGather
        (
            Pstream::linearCommunication(comm),
            Value,
            cop,
            tag,
            comm
        );
    }
    else
    {
        combineGather
        (
            Pstream::treeCommunication(comm),
            Value,
            cop,
            tag,
            comm
        );
    }
}


template <class T>
void Pstream::combineScatter
(
    const List<Pstream::commsStruct>& comms,
    T& Value,
    const int tag,
    const label comm
)
{
    if (Pstream::nProcs(comm) > 1)
    {
        // Get my communication order
        const Pstream::commsStruct& myComm = comms[Pstream::myProcNo(comm)];

        // Receive from up
        if (myComm.above() != -1)
        {
            if (contiguous<T>())
            {
                IPstream::read
                (
                    Pstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<char*>(&Value),
                    sizeof(T),
                    tag,
                    comm
                );
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
                Value = T(fromAbove);
            }

            if (debug > 1)
            {
                Pout<< " received from "
                    << myComm.above() << " data:" << Value << endl;
            }
        }

        // Send to my downstairs neighbours.  Note reverse order (compared to
        // receiving). This is to make sure to send to the critical path
        // (only when using a tree schedule!) first.
        // This is ESI Comms optimisation, v16.06.  HJ, 19/Sep/2016
        forAllReverse (myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            if (debug > 1)
            {
                Pout<< " sending to " << belowID << " data:" << Value << endl;
            }

            if (contiguous<T>())
            {
                OPstream::write
                (
                    Pstream::scheduled,
                    belowID,
                    reinterpret_cast<const char*>(&Value),
                    sizeof(T),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream toBelow(Pstream::scheduled, belowID, 0, tag, comm);
                toBelow << Value;
            }
        }
    }
}


template <class T>
void Pstream::combineScatter
(
    T& Value,
    const int tag,
    const label comm
)
{
    if (Pstream::nProcs(comm) < Pstream::nProcsSimpleSum())
    {
        combineScatter(Pstream::linearCommunication(comm), Value, tag, comm);
    }
    else
    {
        combineScatter(Pstream::treeCommunication(comm), Value, tag, comm);
    }
}


// Same thing but for whole list at a time
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


template <class T, class CombineOp>
void Pstream::listCombineGather
(
    const List<Pstream::commsStruct>& comms,
    List<T>& Values,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (Pstream::nProcs(comm) > 1)
    {
        // Get my communication order
        const commsStruct& myComm = comms[Pstream::myProcNo(comm)];

        // Receive from my downstairs neighbours
        forAll (myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            if (contiguous<T>())
            {
                List<T> receivedValues(Values.size());

                IPstream::read
                (
                    Pstream::scheduled,
                    belowID,
                    reinterpret_cast<char*>(receivedValues.begin()),
                    receivedValues.byteSize(),
                    tag,
                    comm
                );

                if (debug > 1)
                {
                    Pout<< " received from "
                        << belowID << " data:" << receivedValues << endl;
                }

                forAll (Values, i)
                {
                    cop(Values[i], receivedValues[i]);
                }
            }
            else
            {
                IPstream fromBelow(Pstream::scheduled, belowID, 0, tag, comm);
                List<T> receivedValues(fromBelow);

                if (debug > 1)
                {
                    Pout<< " received from "
                        << belowID << " data:" << receivedValues << endl;
                }

                forAll (Values, i)
                {
                    cop(Values[i], receivedValues[i]);
                }
            }
        }

        // Send up Value
        if (myComm.above() != -1)
        {
            if (debug > 1)
            {
                Pout<< " sending to " << myComm.above()
                    << " data:" << Values << endl;
            }

            if (contiguous<T>())
            {
                OPstream::write
                (
                    Pstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<const char*>(Values.begin()),
                    Values.byteSize(),
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
                toAbove << Values;
            }
        }
    }
}


template <class T, class CombineOp>
void Pstream::listCombineGather
(
    List<T>& Values,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (Pstream::nProcs(comm) < Pstream::nProcsSimpleSum())
    {
        listCombineGather
        (
            Pstream::linearCommunication(comm),
            Values,
            cop,
            tag,
            comm
        );
    }
    else
    {
        listCombineGather
        (
            Pstream::treeCommunication(comm),
            Values,
            cop,
            tag,
            comm
        );
    }
}


template <class T>
void Pstream::listCombineScatter
(
    const List<Pstream::commsStruct>& comms,
    List<T>& Values,
    const int tag,
    const label comm
)
{
    if (Pstream::nProcs(comm) > 1)
    {
        // Get my communication order
        const Pstream::commsStruct& myComm = comms[Pstream::myProcNo(comm)];

        // Receive from up
        if (myComm.above() != -1)
        {
            if (contiguous<T>())
            {
                IPstream::read
                (
                    Pstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<char*>(Values.begin()),
                    Values.byteSize(),
                    tag,
                    comm
                );
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
                fromAbove >> Values;
            }

            if (debug > 1)
            {
                Pout<< " received from "
                    << myComm.above() << " data:" << Values << endl;
            }
        }

        // Send to my downstairs neighbours.  Note reverse order (compared to
        // receiving). This is to make sure to send to the critical path
        // (only when using a tree schedule!) first.
        // This is ESI Comms optimisation, v16.06.  HJ, 19/Sep/2016
        forAllReverse (myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            if (debug > 1)
            {
                Pout<< " sending to " << belowID << " data:" << Values << endl;
            }

            if (contiguous<T>())
            {
                OPstream::write
                (
                    Pstream::scheduled,
                    belowID,
                    reinterpret_cast<const char*>(Values.begin()),
                    Values.byteSize(),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream toBelow(Pstream::scheduled, belowID, 0, tag, comm);
                toBelow << Values;
            }
        }
    }
}


template <class T>
void Pstream::listCombineScatter
(
    List<T>& Values,
    const int tag,
    const label comm
)
{
    if (Pstream::nProcs(comm) < Pstream::nProcsSimpleSum())
    {
        listCombineScatter
        (
            Pstream::linearCommunication(),
            Values,
            tag,
            comm
        );
    }
    else
    {
        listCombineScatter
        (
            Pstream::treeCommunication(),
            Values,
            tag,
            comm
        );
    }
}


// Same thing but for sparse list (map)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template <class Container, class CombineOp>
void Pstream::mapCombineGather
(
    const List<Pstream::commsStruct>& comms,
    Container& Values,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (Pstream::nProcs(comm) > 1)
    {
        // Get my communication order
        const commsStruct& myComm = comms[Pstream::myProcNo(comm)];

        // Receive from my downstairs neighbours
        forAll (myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            IPstream fromBelow(Pstream::scheduled, belowID, 0, tag, comm);
            Container receivedValues(fromBelow);

            if (debug > 1)
            {
                Pout<< " received from "
                    << belowID << " data:" << receivedValues << endl;
            }

            for
            (
                typename Container::const_iterator slaveIter =
                    receivedValues.begin();
                slaveIter != receivedValues.end();
                ++slaveIter
            )
            {
                typename Container::iterator
                    masterIter = Values.find(slaveIter.key());

                if (masterIter != Values.end())
                {
                    cop(masterIter(), slaveIter());
                }
                else
                {
                    Values.insert(slaveIter.key(), slaveIter());
                }
            }
        }

        // Send up Value
        if (myComm.above() != -1)
        {
            if (debug > 1)
            {
                Pout<< " sending to " << myComm.above()
                    << " data:" << Values << endl;
            }

            OPstream toAbove(Pstream::scheduled, myComm.above(), 0, tag, comm);
            toAbove << Values;
        }
    }
}


template <class Container, class CombineOp>
void Pstream::mapCombineGather
(
    Container& Values,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (Pstream::nProcs(comm) < Pstream::nProcsSimpleSum())
    {
        mapCombineGather
        (
            Pstream::linearCommunication(),
            Values,
            cop,
            tag,
            comm
        );
    }
    else
    {
        mapCombineGather
        (
            Pstream::treeCommunication(),
            Values,
            cop,
            tag,
            comm
        );
    }
}


template <class Container>
void Pstream::mapCombineScatter
(
    const List<Pstream::commsStruct>& comms,
    Container& Values,
    const int tag,
    const label comm
)
{
    if (Pstream::nProcs(comm) > 1)
    {
        // Get my communication order
        const Pstream::commsStruct& myComm = comms[Pstream::myProcNo(comm)];

        // Receive from up
        if (myComm.above() != -1)
        {
            IPstream fromAbove
            (
                Pstream::scheduled,
                myComm.above(),
                0,
                tag,
                comm
            );
            fromAbove >> Values;

            if (debug > 1)
            {
                Pout<< " received from "
                    << myComm.above() << " data:" << Values << endl;
            }
        }

        // Send to my downstairs neighbours.  Note reverse order (compared to
        // receiving). This is to make sure to send to the critical path
        // (only when using a tree schedule!) first.
        // This is ESI Comms optimisation, v16.06.  HJ, 19/Sep/2016
        forAllReverse (myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            if (debug > 1)
            {
                Pout<< " sending to " << belowID << " data:" << Values << endl;
            }

            OPstream toBelow(Pstream::scheduled, belowID, 0, tag, comm);
            toBelow << Values;
        }
    }
}


template <class Container>
void Pstream::mapCombineScatter
(
    Container& Values,
    const int tag,
    const label comm
)
{
    if (Pstream::nProcs(comm) < Pstream::nProcsSimpleSum())
    {
        mapCombineScatter
        (
            Pstream::linearCommunication(),
            Values,
            tag,
            comm
        );
    }
    else
    {
        mapCombineScatter
        (
            Pstream::treeCommunication(),
            Values,
            tag,
            comm
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
