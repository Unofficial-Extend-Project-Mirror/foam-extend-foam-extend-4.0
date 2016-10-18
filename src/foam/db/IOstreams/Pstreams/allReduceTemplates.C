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

    // Removed send-received loop: use Allreduce instead.
    // HJ, 8/Oct/2016

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::allReduce\n"
            "(\n"
            "    Type& Value,\n"
            "    int MPICount,\n"
            "    MPI_Datatype MPIType,\n"
            "    MPI_Op MPIOp,\n"
            "    const BinaryOp& bop,\n"
            "    const int tag,\n"
            "    const label comm\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

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


// ************************************************************************* //
