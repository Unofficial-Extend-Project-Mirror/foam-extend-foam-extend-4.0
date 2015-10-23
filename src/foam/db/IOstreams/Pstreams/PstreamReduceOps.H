/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#ifndef PstreamReduceOps_H
#define PstreamReduceOps_H

#include "Pstream.H"
#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Reduce operation with user specified communication schedule
template <class T, class BinaryOp>
void reduce
(
    const List<Pstream::commsStruct>& comms,
    T& Value,
    const BinaryOp& bop
)
{
    Pstream::gather(comms, Value, bop);
    Pstream::scatter(comms, Value);
}


// Reduce using either linear or tree communication schedule
template <class T, class BinaryOp>
void reduce
(
    T& Value,
    const BinaryOp& bop
)
{
    if (Pstream::nProcs() < Pstream::nProcsSimpleSum())
    {
        reduce(Pstream::linearCommunication(), Value, bop);
    }
    else
    {
        reduce(Pstream::treeCommunication(), Value, bop);
    }
}


// Reduce using either linear or tree communication schedule
template <class T, class BinaryOp>
T returnReduce
(
    const T& Value,
    const BinaryOp& bop
)
{
    T WorkValue(Value);

    if (Pstream::nProcs() < Pstream::nProcsSimpleSum())
    {
        reduce(Pstream::linearCommunication(), WorkValue, bop);
    }
    else
    {
        reduce(Pstream::treeCommunication(), WorkValue, bop);
    }

    return WorkValue;
}


// Insist there is a specialisation for the reduction of a scalar
void reduce(scalar& Value, const sumOp<scalar>& bop);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
