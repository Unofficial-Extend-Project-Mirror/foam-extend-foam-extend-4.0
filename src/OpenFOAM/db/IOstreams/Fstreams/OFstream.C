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

\*---------------------------------------------------------------------------*/

#include "OFstream.H"
#include "OSspecific.H"
#include "gzstream.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(OFstream, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

OFstreamAllocator::OFstreamAllocator
(
    const fileName& pathname,
    ios_base::openmode mode,
    IOstream::compressionType compression
)
:
    ofPtr_(NULL)
{
    if (!pathname.size())
    {
        if (OFstream::debug)
        {
            Info
                << "OFstreamAllocator::OFstreamAllocator"
                   "(const fileName& pathname) : "
                   "can't open null file "
                << endl;
        }
    }

    if (compression == IOstream::COMPRESSED)
    {
        if (file(pathname))
        {
            rm(pathname);
        }

        ofPtr_ = new ogzstream((pathname + ".gz").c_str());
    }
    else
    {
        if (file(pathname + ".gz"))
        {
            rm(pathname + ".gz");
        }

        ofPtr_ = new ofstream(pathname.c_str());
    }
}


OFstreamAllocator::~OFstreamAllocator()
{
    delete ofPtr_;
}


ostream& OFstreamAllocator::stdStream()
{
    if (!ofPtr_)
    {
        FatalErrorIn("OFstreamAllocator::stdStream()")
            << "No stream allocated." << abort(FatalError);
    }
    return *ofPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

OFstream::OFstream
(
    const fileName& pathname,
    ios_base::openmode mode,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    OFstreamAllocator(pathname, mode, compression),
    OSstream(*ofPtr_, "OFstream.sinkFile_", format, version, compression),
    pathname_(pathname)
{
    setClosed();

    setState(ofPtr_->rdstate());

    if (!good())
    {
        if (debug)
        {
            Info<< "IFstream::IFstream(const fileName& pathname,"
                   "streamFormat format=ASCII,"
                   "versionNumber version=currentVersion) : "
                   "couldn't open File for input\n"
                   "in stream " << info() << Foam::endl;
        }

        setBad();
    }
    else
    {
        setOpened();
    }

    lineNumber_ = 1;
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

OFstream::~OFstream()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void OFstream::print(Ostream& os) const
{
    // Print File data
    os  << "    OFstream: ";
    OSstream::print(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
