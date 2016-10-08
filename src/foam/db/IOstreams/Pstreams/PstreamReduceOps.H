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

#ifndef PstreamReduceOps_H
#define PstreamReduceOps_H

#include "Pstream.H"
#include "ops.H"
#include "vector2D.H"

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
    const BinaryOp& bop,
    const int tag,
    const label comm
)
{
    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << endl;
        error::printStack(Pout);
    }

    Pstream::gather(comms, Value, bop, tag, comm);
    Pstream::scatter(comms, Value, tag, comm);
}


// Reduce using either linear or tree communication schedule
template <class T, class BinaryOp>
void reduce
(
    T& Value,
    const BinaryOp& bop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
)
{
    if (Pstream::nProcs(comm) < Pstream::nProcsSimpleSum())
    {
        reduce(Pstream::linearCommunication(comm), Value, bop, tag, comm);
    }
    else
    {
        reduce(Pstream::treeCommunication(comm), Value, bop, tag, comm);
    }
}


// Reduce using either linear or tree communication schedule
template <class T, class BinaryOp>
T returnReduce
(
    const T& Value,
    const BinaryOp& bop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
)
{
    T WorkValue(Value);

    if (Pstream::nProcs(comm) < Pstream::nProcsSimpleSum())
    {
        reduce(Pstream::linearCommunication(comm), WorkValue, bop, tag, comm);
    }
    else
    {
        reduce(Pstream::treeCommunication(comm), WorkValue, bop, tag, comm);
    }

    return WorkValue;
}


// Reduce with sum of both value and count (for averaging)
template<class T>
void sumReduce
(
    T& Value,
    label& Count,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
)
{
    reduce(Value, sumOp<T>(), tag, comm);
    reduce(Count, sumOp<label>(), tag, comm);
}


// Non-blocking version of reduce. Sets request.
template<class T, class BinaryOp>
void reduce
(
    T& Value,
    const BinaryOp& bop,
    const int tag,
    const label comm,
    label& request
)
{
    notImplemented
    (
        "reduce(T&, const BinaryOp&, const int, const label, label&)"
    );
}


// Insist there are specialisations for the common reductions of bool
// See implementation notes on casting of a C++ bool into MPI_INT in
// Pstream.C  HJ, 23/Sep/2016

void reduce
(
    bool& Value,
    const andOp<bool>& bop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
);

void reduce
(
    bool& Value,
    const orOp<bool>& bop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
);


// Insist there are specialisations for the common reductions of labels

void reduce
(
    label& Value,
    const sumOp<label>& bop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
);

void reduce
(
    label& Value,
    const minOp<label>& bop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
);

void reduce
(
    label& Value,
    const maxOp<label>& bop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
);

// Insist there are specialisations for the common reductions of scalars

void reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
);

void reduce
(
    scalar& Value,
    const minOp<scalar>& bop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
);

void reduce
(
    scalar& Value,
    const maxOp<scalar>& bop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
);

void reduce
(
    vector2D& Value,
    const sumOp<vector2D>& bop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
);

void sumReduce
(
    scalar& Value,
    label& Count,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
);

void reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label comm,
    label& request
);


// Insist there are specialisations for the common reductions of
// lists of labels.  Note: template function specialisation must be the
// exact match on argument types.  HJ, 8/Oct/2016

void reduce
(
    List<label>& Value,
    const sumOp<List<label> >& bop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
);

void reduce
(
    List<label>& Value,
    const minOp<List<label> >& bop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
);

void reduce
(
    List<label>& Value,
    const maxOp<List<label> >& bop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
