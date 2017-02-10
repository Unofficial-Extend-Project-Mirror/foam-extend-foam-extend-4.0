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
    Read token and binary block from IPstream

\*---------------------------------------------------------------------------*/

#include "mpi.h"

#include "IPstream.H"
#include "PstreamGlobals.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::IPstream::IPstream
(
    const commsTypes commsType,
    const int fromProcNo,
    const label bufSize,
    const int tag,
    const label comm,
    streamFormat format,
    versionNumber version
)
:
    Pstream(commsType, bufSize),
    Istream(format, version),
    fromProcNo_(fromProcNo),
    tag_(tag),
    comm_(comm),
    messageSize_(0)
{
    setOpened();
    setGood();

    MPI_Status status;

    // If the buffer size is not specified, probe the incoming message
    // and set it
    if (!bufSize)
    {
        MPI_Probe
        (
            fromProcNo_,
            tag_,
            PstreamGlobals::MPICommunicators_[comm_],
            &status
        );
        MPI_Get_count(&status, MPI_BYTE, &messageSize_);

        buf_.setSize(messageSize_);
    }

    messageSize_ = IPstream::read
    (
        commsType,
        fromProcNo_,
        buf_.begin(),
        buf_.size(),
        tag_,
        comm_
    );

    if (!messageSize_)
    {
        FatalErrorIn
        (
            "IPstream::IPstream(const int fromProcNo, "
            "const label bufSize, streamFormat format, versionNumber version)"
        )   << "read failed"
            << Foam::abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::IPstream::read
(
    const commsTypes commsType,
    const int fromProcNo,
    char* buf,
    const std::streamsize bufSize,
    const int tag,
    const label comm
)
{
    if (debug)
    {
        Pout<< "IPstream::read : starting read from:" << fromProcNo
            << " tag:" << tag << " comm:" << comm
            << " wanted size:" << label(bufSize)
            << " commsType:" << Pstream::commsTypeNames[commsType]
            << Foam::endl;
    }

    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "IPstream::read : starting read from:" << fromProcNo
            << " tag:" << tag << " comm:" << comm
            << " wanted size:" << label(bufSize)
            << " commsType:" << Pstream::commsTypeNames[commsType]
            << " warnComm:" << Pstream::warnComm
            << Foam::endl;
        error::printStack(Pout);
    }

    if (commsType == blocking || commsType == scheduled)
    {
        MPI_Status status;

        if
        (
            MPI_Recv
            (
                buf,
                bufSize,
                MPI_BYTE,
                fromProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[comm],
                &status
            )
        )
        {
            FatalErrorIn
            (
                "IPstream::read"
                "(const int fromProcNo, char* buf, std::streamsize bufSize)"
            )   << "MPI_Recv cannot receive incomming message"
                << Foam::abort(FatalError);

            return 0;
        }


        // Check size of message read

        label messageSize;
        MPI_Get_count(&status, MPI_BYTE, &messageSize);

        if (messageSize > bufSize)
        {
            FatalErrorIn
            (
                "IPstream::read"
                "(const int fromProcNo, char* buf, std::streamsize bufSize)"
            )   << "buffer (" << label(bufSize)
                << ") not large enough for incomming message ("
                << messageSize << ')'
                << Foam::abort(FatalError);
        }

        return messageSize;
    }
    else if (commsType == nonBlocking)
    {
        MPI_Request request;

        if
        (
            MPI_Irecv
            (
                buf,
                bufSize,
                MPI_BYTE,
                fromProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[comm],
                &request
            )
        )
        {
            FatalErrorIn
            (
                "IPstream::read"
                "(const int fromProcNo, char* buf, std::streamsize bufSize)"
            )   << "MPI_Recv cannot start non-blocking receive"
                << Foam::abort(FatalError);

            return 0;
        }

        if (debug)
        {
            Pout<< "IPstream::read : started read from:" << fromProcNo
                << " tag:" << tag << " read size:" << label(bufSize)
                << " commsType:" << Pstream::commsTypeNames[commsType]
                << " request:" << PstreamGlobals::outstandingRequests_.size()
                << Foam::endl;
        }

        PstreamGlobals::outstandingRequests_.append(request);

        // Assume the message is completely received.
        return 1;
    }
    else
    {
        FatalErrorIn
        (
            "IPstream::read"
            "(const int fromProcNo, char* buf, std::streamsize bufSize)"
        )   << "Unsupported communications type " << commsType
            << Foam::abort(FatalError);

        return 0;
    }
}


// ************************************************************************* //
