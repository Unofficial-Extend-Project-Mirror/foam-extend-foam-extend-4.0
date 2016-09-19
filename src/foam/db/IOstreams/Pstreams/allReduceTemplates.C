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

#include "allReduce.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type, class BinaryOp>
void Foam::allReduce
(
    Type& Value,
    int MPICount,
    MPI_Datatype MPIType,
    MPI_Op MPIOp,
    const BinaryOp& bop,
    const int tag,
    const label comm
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (Pstream::nProcs(comm) <= Pstream::nProcsSimpleSum)
    {
        if (Pstream::master(comm))
        {
            for
            (
                int slave = Pstream::firstSlave();
                slave <= Pstream::lastSlave(comm);
                slave++
            )
            {
                Type value;

                if
                (
                    MPI_Recv
                    (
                        &value,
                        MPICount,
                        MPIType,
                        slave,
                        tag,
                        PstreamGlobals::MPICommunicators_[comm],
                        MPI_STATUS_IGNORE
                    )
                )
                {
                    FatalErrorIn
                    (
                        "void Foam::allReduce\n"
                        "(\n"
                        "    Type&,\n"
                        "    int,\n"
                        "    MPI_Datatype,\n"
                        "    MPI_Op,\n"
                        "    const BinaryOp&,\n"
                        "    const int\n"
                        ")\n"
                    )   << "MPI_Recv failed"
                        << Foam::abort(FatalError);
                }

                Value = bop(Value, value);
            }
        }
        else
        {
            if
            (
                MPI_Send
                (
                    &Value,
                    MPICount,
                    MPIType,
                    Pstream::masterNo(),
                    tag,
                    PstreamGlobals::MPICommunicators_[comm]
                )
            )
            {
                FatalErrorIn
                (
                    "void Foam::allReduce\n"
                    "(\n"
                    "    Type&,\n"
                    "    int,\n"
                    "    MPI_Datatype,\n"
                    "    MPI_Op,\n"
                    "    const BinaryOp&,\n"
                    "    const int\n"
                    ")\n"
                )   << "MPI_Send failed"
                    << Foam::abort(FatalError);
            }
        }


        if (Pstream::master(comm))
        {
            for
            (
                int slave = Pstream::firstSlave();
                slave <= Pstream::lastSlave(comm);
                slave++
            )
            {
                if
                (
                    MPI_Send
                    (
                        &Value,
                        MPICount,
                        MPIType,
                        slave,
                        tag,
                        PstreamGlobals::MPICommunicators_[comm]
                    )
                )
                {
                    FatalErrorIn
                    (
                        "void Foam::allReduce\n"
                        "(\n"
                        "    Type&,\n"
                        "    int,\n"
                        "    MPI_Datatype,\n"
                        "    MPI_Op,\n"
                        "    const BinaryOp&,\n"
                        "    const int\n"
                        ")\n"
                    )   << "MPI_Send failed"
                        << Foam::abort(FatalError);
                }
            }
        }
        else
        {
            if
            (
                MPI_Recv
                (
                    &Value,
                    MPICount,
                    MPIType,
                    Pstream::masterNo(),
                    tag,
                    PstreamGlobals::MPICommunicators_[comm],
                    MPI_STATUS_IGNORE
                )
            )
            {
                FatalErrorIn
                (
                    "void Foam::allReduce\n"
                    "(\n"
                    "    Type&,\n"
                    "    int,\n"
                    "    MPI_Datatype,\n"
                    "    MPI_Op,\n"
                    "    const BinaryOp&,\n"
                    "    const int\n"
                    ")\n"
                )   << "MPI_Recv failed"
                    << Foam::abort(FatalError);
            }
        }
    }
    else
    {
        Type sum;

        MPI_Allreduce
        (
            &Value,
            &sum,
            MPICount,
            MPIType,
            MPIOp,
            PstreamGlobals::MPICommunicators_[comm]
        );

        Value = sum;
    }
}


// ************************************************************************* //
