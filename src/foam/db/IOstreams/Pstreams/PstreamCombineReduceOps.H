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

InClass
    Foam::Pstream

Description
    Combination-Reduction operation for a parallel run.  The
    information from all nodes is collected on the master node,
    combined using the given combination function and the result is
    broadcast to all nodes

\*---------------------------------------------------------------------------*/

#ifndef PstreamCombineReduceOps_H
#define PstreamCombineReduceOps_H

#include "Pstream.H"
#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <class T, class CombineOp>
void combineReduce
(
    const List<Pstream::commsStruct>& comms,
    T& Value,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    Pstream::combineGather(comms, Value, cop, tag, comm);
    Pstream::combineScatter(comms, Value, tag, comm);
}


template <class T, class CombineOp>
void combineReduce
(
    T& Value,
    const CombineOp& cop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
)
{
    if (Pstream::nProcs() < Pstream::nProcsSimpleSum())
    {
        Pstream::combineGather
        (
            Pstream::linearCommunication(comm),
            Value,
            cop,
            tag,
            comm
        );

        Pstream::combineScatter
        (
            Pstream::linearCommunication(comm),
            Value,
            tag,
            comm
        );
    }
    else
    {
        Pstream::combineGather
        (
            Pstream::treeCommunication(comm),
            Value,
            cop,
            tag,
            comm
        );

        Pstream::combineScatter
        (
            Pstream::treeCommunication(comm),
            Value,
            tag,
            comm
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
