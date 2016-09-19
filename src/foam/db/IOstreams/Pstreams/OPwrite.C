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
    Write primitive and binary block from OPstream

\*---------------------------------------------------------------------------*/

#include "mpi.h"

#include "OPstream.H"
#include "PstreamGlobals.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::OPstream::~OPstream()
{
    if
    (
       !write
        (
            commsType_,
            toProcNo_,
            buf_.begin(),
            bufPosition_,
            tag_,
            comm_
        )
    )
    {
        FatalErrorIn("OPstream::~OPstream()")
            << "MPI_Bsend cannot send outgoing message"
            << Foam::abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::OPstream::write
(
    const commsTypes commsType,
    const int toProcNo,
    const char* buf,
    const std::streamsize bufSize,
    const int tag,
    const label comm
)
{
    if (debug)
    {
        Pout<< "OPstream::write : starting write to:" << toProcNo
            << " tag:" << tag
            << " comm:" << comm << " size:" << label(bufSize)
            << " commsType:" << Pstream::commsTypeNames[commsType]
            << Foam::endl;
    }
    if (Pstream::warnComm != -1 && comm != Pstream::warnComm)
    {
        Pout<< "OPstream::write : starting write to:" << toProcNo
            << " tag:" << tag
            << " comm:" << comm << " size:" << label(bufSize)
            << " commsType:" << Pstream::commsTypeNames[commsType]
            << " warnComm:" << Pstream::warnComm
            << Foam::endl;

        error::printStack(Pout);
    }


    PstreamGlobals::checkCommunicator(comm, toProcNo);

    bool transferFailed = true;

    if (commsType == blocking)
    {
        transferFailed = MPI_Bsend
        (
            const_cast<char*>(buf),
            bufSize,
            MPI_BYTE,
            toProcNo,   //procID(toProcNo),
            tag,
            PstreamGlobals::MPICommunicators_[comm] // MPI_COMM_WORLD
        );

        if (debug)
        {
            Pout<< "OPstream::write : finished write to:" << toProcNo
                << " tag:" << tag << " size:" << label(bufSize)
                << " commsType:" << Pstream::commsTypeNames[commsType]
                << Foam::endl;
        }
    }
    else if (commsType == scheduled)
    {
        transferFailed = MPI_Send
        (
            const_cast<char*>(buf),
            bufSize,
            MPI_BYTE,
            toProcNo,   //procID(toProcNo),
            tag,
            PstreamGlobals::MPICommunicators_[comm] // MPI_COMM_WORLD
        );

        if (debug)
        {
            Pout<< "OPstream::write : finished write to:" << toProcNo
                << " tag:" << tag << " size:" << label(bufSize)
                << " commsType:" << Pstream::commsTypeNames[commsType]
                << Foam::endl;
        }
    }
    else if (commsType == nonBlocking)
    {
        MPI_Request request;

        transferFailed = MPI_Isend
        (
            const_cast<char*>(buf),
            bufSize,
            MPI_BYTE,
            toProcNo,   //procID(toProcNo),
            tag,
            PstreamGlobals::MPICommunicators_[comm],// MPI_COMM_WORLD,
            &request
        );

        if (debug)
        {
            Pout<< "OPstream::write : started write to:" << toProcNo
                << " tag:" << tag << " size:" << label(bufSize)
                << " commsType:" << Pstream::commsTypeNames[commsType]
                << " request:" << PstreamGlobals::outstandingRequests_.size()
                << Foam::endl;
        }

        PstreamGlobals::outstandingRequests_.append(request);
    }
    else
    {
        FatalErrorIn
        (
            "OPstream::write"
            "(const int fromProcNo, char* buf, std::streamsize bufSize)"
        )   << "Unsupported communications type " << commsType
            << Foam::abort(FatalError);
    }

    return !transferFailed;
}


// ************************************************************************* //
