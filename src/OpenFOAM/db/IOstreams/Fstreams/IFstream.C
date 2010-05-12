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

#include "IFstream.H"
#include "OSspecific.H"
#include "gzstream.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(IFstream, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

IFstreamAllocator::IFstreamAllocator(const fileName& pathname)
:
    ifPtr_(NULL),
    compression_(IOstream::UNCOMPRESSED)
{
    if (!pathname.size())
    {
        if (IFstream::debug)
        {
            Info<< "IFstreamAllocator::IFstreamAllocator"
                   "(const fileName& pathname) : "
                   "can't open null file "
                << endl;
        }
    }

    ifPtr_ = new ifstream(pathname.c_str());

    // If the file is compressed, decompress it before reading.
    if (!ifPtr_->good() && file(pathname + ".gz"))
    {
        if (IFstream::debug)
        {
            Info<< "IFstreamAllocator::IFstreamAllocator"
                   "(const fileName& pathname) : "
                   "decompressing " << pathname + ".gz"
                << endl;
        }

        delete ifPtr_;

        ifPtr_ = new igzstream((pathname + ".gz").c_str());

        if (ifPtr_->good())
        {
            compression_ = IOstream::COMPRESSED;
        }
    }
}


IFstreamAllocator::~IFstreamAllocator()
{
    delete ifPtr_;
}


istream& IFstreamAllocator::stdStream()
{
    if (!ifPtr_)
    {
        FatalErrorIn("IFstreamAllocator::stdStream()")
            << "No stream allocated." << abort(FatalError);
    }
    return *ifPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

IFstream::IFstream
(
    const fileName& pathname,
    streamFormat format,
    versionNumber version
)
:
    IFstreamAllocator(pathname),
    ISstream
    (
        *ifPtr_,
        "IFstream.sourceFile_",
        format,
        version,
        IFstreamAllocator::compression_
    ),
    pathname_(pathname)
{
    setClosed();

    setState(ifPtr_->rdstate());
                
    if (!good())
    {
        if (debug)
        {
            Info<< "IFstream::IFstream(const fileName& pathname,"
                   "streamFormat format=ASCII,"
                   "versionNumber version=currentVersion) : "
                   "couldn't open File for input"
                << endl << info() << endl;
        }

        setBad();
    }
    else
    {
        setOpened();
    }
    
    lineNumber_ = 1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

IFstream::~IFstream()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void IFstream::print(Ostream& os) const
{
    // Print File data
    os  << "IFstream: ";
    ISstream::print(os);
}


//- Return a non-const reference to const Istream
//  Needed for read-constructors where the stream argument is temporary:
//  e.g. thing thisThing(IFstream("thingFileName")());
IFstream& IFstream::operator()() const
{
    if (!good())
    {
        if (!file(pathname_) && !file(pathname_ + ".gz"))
        {
            FatalIOErrorIn("IFstream::operator()", *this)
                << "file " << pathname_ << " does not exist"
                << exit(FatalIOError);
        }
        else
        {
            check("IFstream::operator()");
            FatalIOError.exit();
        }
    }

    return const_cast<IFstream&>(*this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
