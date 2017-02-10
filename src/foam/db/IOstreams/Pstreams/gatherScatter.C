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
    The gathered data will be a single value constructed from the values
    on individual processors using a user-specified operator.

\*---------------------------------------------------------------------------*/

#include "OPstream.H"
#include "IPstream.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T, class BinaryOp>
void Pstream::gather
(
    const List<Pstream::commsStruct>& comms,
    T& Value,
    const BinaryOp& bop,
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
        // Get my communication order
        const commsStruct& myComm = comms[Pstream::myProcNo(comm)];

        // Receive from my downstairs neighbours
        forAll (myComm.below(), belowI)
        {
            T value;

            if (contiguous<T>())
            {
                IPstream::read
                (
                    Pstream::scheduled,
                    myComm.below()[belowI],
                    reinterpret_cast<char*>(&value),
                    sizeof(T),
                    tag,
                    comm
                );
            }
            else
            {
                IPstream fromBelow
                (
                    Pstream::scheduled,
                    myComm.below()[belowI],
                    0,
                    tag,
                    comm
                );
                fromBelow >> value;
            }

            Value = bop(Value, value);
        }

        // Send up Value
        if (myComm.above() != -1)
        {
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


template <class T, class BinaryOp>
void Pstream::gather
(
    T& Value,
    const BinaryOp& bop,
    const int tag,
    const label comm
)
{
    if (Pstream::nProcs(comm) < Pstream::nProcsSimpleSum())
    {
        gather(Pstream::linearCommunication(comm), Value, bop, tag, comm);
    }
    else
    {
        gather(Pstream::treeCommunication(comm), Value, bop, tag, comm);
    }
}


template <class T>
void Pstream::scatter
(
    const List<Pstream::commsStruct>& comms,
    T& Value,
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
        // Get my communication order
        const commsStruct& myComm = comms[Pstream::myProcNo(comm)];

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
                fromAbove >> Value;
            }
        }

        // Send to my downstairs neighbours.  Note reverse order (compared to
        // receiving). This is to make sure to send to the critical path
        // (only when using a tree schedule!) first.
        // This is ESI Comms optimisation, v16.06.  HJ, 19/Sep/2016
        forAllReverse (myComm.below(), belowI)
        {
            if (contiguous<T>())
            {
                OPstream::write
                (
                    Pstream::scheduled,
                    myComm.below()[belowI],
                    reinterpret_cast<const char*>(&Value),
                    sizeof(T),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream toBelow
                (
                    Pstream::scheduled,
                    myComm.below()[belowI],
                    0,
                    tag,
                    comm
                );
                toBelow << Value;
            }
        }
    }
}


template <class T>
void Pstream::scatter(T& Value, const int tag, const label comm)
{
    if (Pstream::nProcs(comm) < Pstream::nProcsSimpleSum())
    {
        scatter(Pstream::linearCommunication(comm), Value, tag, comm);
    }
    else
    {
        scatter(Pstream::treeCommunication(comm), Value, tag, comm);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
