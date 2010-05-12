/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Read token and binary block from IPstream using pvm.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "IPstream.H"

#include <pvm3.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

IPstream::IPstream
(
    const int fromProcNo,
    const label bufSize,
    streamFormat format,
    versionNumber version
)
:
    Pstream(bufSize),
    Istream(format, version),
    fromProcNo_(fromProcNo),
    messageSize_(0)
{
    setOpened();
    setGood();

    int bufid, tag, tid;

    // If the buffer size is not specified then probe the incomming message

    if (!bufSize)
    {
        // Probe read buffer until message arrives.
        while (!(bufid = pvm_probe(procID(fromProcNo_), msgType())));

        // When the message arrives find its size
        pvm_bufinfo(bufid, &messageSize_, &tag, &tid);

        // Resize buffer to message size
        buf_.setSize(messageSize_);
    }


    // Read message into buffer

    if
    (
        pvm_precv
        (
            procID(fromProcNo_),
            msgType(),
            buf_.begin(),
            buf_.size(),
            PVM_BYTE,
            &tid, &tag, &messageSize_
        ) != PvmOk
    )
    {
        FatalErrorIn("IPstream::IPstream(const int fromProcNo)")
            << "pvm_precv cannot receive incomming message"
            << ::abort;
    }


    // Check size of message read

    if (messageSize_ > buf_.size())
    {
        FatalErrorIn("IPstream::IPstream(const int fromProcNo)")
            << "buffer (" << buf_.size()
            << ") not large enough for incomming message ("
            << messageSize_ << ')'
            << ::abort;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
