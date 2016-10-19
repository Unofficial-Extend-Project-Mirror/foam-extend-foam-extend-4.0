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

#include "mpi.h"

#include "label.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "allReduce.H"

// Check type of label for use in MPI calls
#if defined(WM_INT)
#   define MPI_LABEL MPI_INT
#elif defined(WM_LONG)
#   define MPI_LABEL MPI_LONG
#elif defined(WM_LLONG)
#   define MPI_LABEL MPI_LONG_LONG
#endif

// Check type of scalar for use in MPI calls
#if defined(WM_SP)
#   define MPI_SCALAR MPI_FLOAT
#elif defined(WM_DP)
#   define MPI_SCALAR MPI_DOUBLE
#elif defined(WM_LDP)
#   define MPI_SCALAR MPI_LONG_DOUBLE
#endif


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Optimizing template specialisations using MPI_REDUCE
void Foam::reduce
(
    bool& Value,
    const andOp<bool>& bop,
    const int tag,
    const label comm
)
{
    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << " warnComm:" << Pstream::warnComm
            << endl;
        error::printStack(Pout);
    }

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::reduce\n"
            "(\n"
            "    bool& Value,\n"
            "    const andOp<bool>& bop,\n"
            "    const int tag,\n"
            "    const label comm\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

    // Note: C++ bool is a type separate from C and cannot be cast
    // For safety and compatibility with compilers, convert bool to int
    // to comply with MPI types.  HJ, 23/Sep/2016

    int intBool = 0;

    if (Value)
    {
        intBool = 1;
    }

    allReduce(intBool, 1, MPI_INT, MPI_LAND, bop, tag, comm);

    if (intBool > 0)
    {
        Value = true;
    }
    else
    {
        Value = false;
    }
}


void Foam::reduce
(
    bool& Value,
    const orOp<bool>& bop,
    const int tag,
    const label comm
)
{
    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << " warnComm:" << Pstream::warnComm
            << endl;
        error::printStack(Pout);
    }

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::reduce\n"
            "(\n"
            "    bool& Value,\n"
            "    const orOp<bool>& bop,\n"
            "    const int tag,\n"
            "    const label comm\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

    // Note: C++ bool is a type separate from C and cannot be cast
    // For safety and compatibility with compilers, convert bool to int
    // to comply with MPI types.  HJ, 23/Sep/2016

    int intBool = 0;

    if (Value)
    {
        intBool = 1;
    }

    allReduce(intBool, 1, MPI_INT, MPI_LOR, bop, tag, comm);

    if (intBool > 0)
    {
        Value = true;
    }
    else
    {
        Value = false;
    }
}


void Foam::reduce
(
    label& Value,
    const minOp<label>& bop,
    const int tag,
    const label comm
)
{
    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << " warnComm:" << Pstream::warnComm
            << endl;
        error::printStack(Pout);
    }

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::reduce\n"
            "(\n"
            "    label& Value,\n"
            "    const minOp<label>& bop,\n"
            "    const int tag,\n"
            "    const label comm\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

    allReduce(Value, 1, MPI_LABEL, MPI_MIN, bop, tag, comm);
}


void Foam::reduce
(
    label& Value,
    const maxOp<label>& bop,
    const int tag,
    const label comm
)
{
    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << " warnComm:" << Pstream::warnComm
            << endl;
        error::printStack(Pout);
    }

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::reduce\n"
            "(\n"
            "    label& Value,\n"
            "    const maxOp<label>& bop,\n"
            "    const int tag,\n"
            "    const label comm\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

    allReduce(Value, 1, MPI_LABEL, MPI_MAX, bop, tag, comm);
}


void Foam::reduce
(
    label& Value,
    const sumOp<label>& bop,
    const int tag,
    const label comm
)
{
    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << " warnComm:" << Pstream::warnComm
            << endl;
        error::printStack(Pout);
    }

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::reduce\n"
            "(\n"
            "    label& Value,\n"
            "    const sumOp<label>& bop,\n"
            "    const int tag,\n"
            "    const label comm\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

    allReduce(Value, 1, MPI_LABEL, MPI_SUM, bop, tag, comm);
}


void Foam::reduce
(
    scalar& Value,
    const minOp<scalar>& bop,
    const int tag,
    const label comm
)
{
    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << " warnComm:" << Pstream::warnComm
            << endl;
        error::printStack(Pout);
    }

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::reduce\n"
            "(\n"
            "    scalar& Value,\n"
            "    const minOp<scalar>& bop,\n"
            "    const int tag,\n"
            "    const label comm\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

    allReduce(Value, 1, MPI_SCALAR, MPI_MIN, bop, tag, comm);
}


void Foam::reduce
(
    scalar& Value,
    const maxOp<scalar>& bop,
    const int tag,
    const label comm
)
{
    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << " warnComm:" << Pstream::warnComm
            << endl;
        error::printStack(Pout);
    }

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::reduce\n"
            "(\n"
            "    scalar& Value,\n"
            "    const maxOp<scalar>& bop,\n"
            "    const int tag,\n"
            "    const label comm\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

    allReduce(Value, 1, MPI_SCALAR, MPI_MAX, bop, tag, comm);
}


void Foam::reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label comm
)
{
    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << " warnComm:" << Pstream::warnComm
            << endl;
        error::printStack(Pout);
    }

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::reduce\n"
            "(\n"
            "    scalar& Value,\n"
            "    const sumOp<scalar>& bop,\n"
            "    const int tag,\n"
            "    const label comm\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

    allReduce(Value, 1, MPI_SCALAR, MPI_SUM, bop, tag, comm);
}


void Foam::reduce
(
    List<label>& Value,
    const minOp<List<label> >& bop,
    const int tag,
    const label comm
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << " warnComm:" << Pstream::warnComm
            << endl;
        error::printStack(Pout);
    }

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::reduce\n"
            "(\n"
            "    List<label>& Value,\n"
            "    const minOp<List<label> >& bop,\n"
            "    const int tag,\n"
            "    const label comm\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

    // Make a copy of the Value in the send buffer so that Value can be
    // used for receive.  HJ, 8/Oct/2016
    labelList send(Value);

    int MPISize = Value.size();

    MPI_Allreduce
    (
        send.begin(),
        Value.begin(),
        MPISize,
        MPI_LABEL,
        MPI_MIN,
        PstreamGlobals::MPICommunicators_[comm]
    );
}


void Foam::reduce
(
    List<label>& Value,
    const maxOp<List<label> >& bop,
    const int tag,
    const label comm
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << " warnComm:" << Pstream::warnComm
            << endl;
        error::printStack(Pout);
    }

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::reduce\n"
            "(\n"
            "    List<label>& Value,\n"
            "    const maxOp<List<label> >& bop,\n"
            "    const int tag,\n"
            "    const label comm\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

    // Make a copy of the Value in the send buffer so that Value can be
    // used for receive.  HJ, 8/Oct/2016
    labelList send(Value);

    int MPISize = Value.size();

    MPI_Allreduce
    (
        send.begin(),
        Value.begin(),
        MPISize,
        MPI_LABEL,
        MPI_MAX,
        PstreamGlobals::MPICommunicators_[comm]
    );
}


void Foam::reduce
(
    List<label>& Value,
    const sumOp<List<label> >& bop,
    const int tag,
    const label comm
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << " warnComm:" << Pstream::warnComm
            << endl;
        error::printStack(Pout);
    }

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::reduce\n"
            "(\n"
            "    List<label>& Value,\n"
            "    const sumOp<List<label> >& bop,\n"
            "    const int tag,\n"
            "    const label comm\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

    // Make a copy of the Value in the send buffer so that Value can be
    // used for receive.  HJ, 8/Oct/2016
    labelList send(Value);

    int MPISize = Value.size();

    MPI_Allreduce
    (
        send.begin(),
        Value.begin(),
        MPISize,
        MPI_LABEL,
        MPI_SUM,
        PstreamGlobals::MPICommunicators_[comm]
    );
}


void Foam::reduce
(
    vector2D& Value,
    const sumOp<vector2D>& bop,
    const int tag,
    const label comm
)
{
    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << " warnComm:" << Pstream::warnComm
            << endl;
        error::printStack(Pout);
    }

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::reduce\n"
            "(\n"
            "    vector2D& Value,\n"
            "    const sumOp<vector2D>& bop,\n"
            "    const int tag,\n"
            "    const label comm\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

    allReduce(Value, 2, MPI_SCALAR, MPI_SUM, bop, tag, comm);
}


void Foam::sumReduce
(
    scalar& Value,
    label& Count,
    const int tag,
    const label comm
)
{
    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << " warnComm:" << Pstream::warnComm
            << endl;
        error::printStack(Pout);
    }

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::sumReduce\n"
            "(\n"
            "    scalar& Value,\n"
            "    label& Count,\n"
            "    const int tag,\n"
            "    const label comm\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

    vector2D twoScalars(Value, scalar(Count));
    reduce(twoScalars, sumOp<vector2D>(), tag, comm);

    Value = twoScalars.x();
    Count = twoScalars.y();
}


void Foam::reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label comm,
    label& requestID
)
{
    if (!Pstream::parRun())
    {
        return;
    }

#   ifdef FULLDEBUG
    // Check for processors that are not in the communicator
    if (Pstream::myProcNo(comm) == -1)
    {
        FatalErrorIn
        (
            "void Foam::reduce\n"
            "(\n"
            "    scalar& Value,\n"
            "    const sumOp<scalar>& bop,\n"
            "    const int tag,\n"
            "    const label comm,\n"
            "    label& requestID\n"
            ")"
        )   << "Reduce called on the processor which is not a member "
            << "of comm.  This is not allowed"
            << abort(FatalError);
    }
#   endif

#ifdef MPIX_COMM_TYPE_SHARED
    // Assume mpich2 with non-blocking collectives extensions. Once mpi3
    // is available this will change.
    MPI_Request request;
    scalar v = Value;
    MPIX_Ireduce
    (
        &v,
        &Value,
        1,
        MPI_SCALAR,
        MPI_SUM,
        0,              // root
        PstreamGlobals::MPICommunicators_[comm],
        &request
    );

    requestID = PstreamGlobals::outstandingRequests_.size();
    PstreamGlobals::outstandingRequests_.append(request);

    if (debug)
    {
        Pout<< "Pstream::allocateRequest for non-blocking reduce"
            << " : request:" << requestID
            << endl;
    }
#else
    // Non-blocking not yet implemented in mpi
    reduce(Value, bop, tag, comm);

    requestID = -1;
#endif
}


// ************************************************************************* //
