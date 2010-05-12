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

#include "IOobject.H"
#include "Time.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IOobject, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOobject::IOobject
(
    const word& name,
    const fileName& instance,
    const objectRegistry& registry,
    readOption ro,
    writeOption wo,
    bool registerObject
)
:
    name_(name),
    headerClassName_(typeName),
    note_(),
    instance_(instance),
    local_(),
    db_(registry),
    rOpt_(ro),
    wOpt_(wo),
    registerObject_(registerObject),
    objState_(GOOD)
{
    if (objectRegistry::debug)
    {
        Info<< "Constructing IOobject called " << name_
            << " of type " << headerClassName_
            << endl;
    }
}


Foam::IOobject::IOobject
(
    const word& name,
    const fileName& instance,
    const fileName& local,
    const objectRegistry& registry,
    readOption ro,
    writeOption wo,
    bool registerObject
)
:
    name_(name),
    headerClassName_(typeName),
    note_(),
    instance_(instance),
    local_(local),
    db_(registry),
    rOpt_(ro),
    wOpt_(wo),
    registerObject_(registerObject),
    objState_(GOOD)
{
    if (objectRegistry::debug)
    {
        Info<< "Constructing IOobject called " << name_
            << " of type " << headerClassName_
            << endl;
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::IOobject::~IOobject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::objectRegistry& Foam::IOobject::db() const
{
    return db_;
}


const Foam::Time& Foam::IOobject::time() const
{
    return db_.time();
}


const Foam::fileName& Foam::IOobject::caseName() const
{
    return time().caseName();
}


const Foam::fileName& Foam::IOobject::rootPath() const
{
    return time().rootPath();
}


Foam::fileName Foam::IOobject::path() const
{
    return rootPath()/caseName()/instance()/db_.dbDir()/local();
}


Foam::fileName Foam::IOobject::path
(
    const word& instance,
    const fileName& local
) const
{
    return rootPath()/caseName()/instance/db_.dbDir()/local;
}


Foam::fileName Foam::IOobject::filePath() const
{
    fileName path = this->path();
    fileName objectPath = path/name();

    if (file(objectPath))
    {
        return objectPath;
    }
    else
    {
        if
        (
            time().processorCase()
         && (
                instance() == time().system()
             || instance() == time().constant()
            )
        )
        {
            fileName parentObjectPath =
                rootPath()/caseName()
               /".."/instance()/db_.dbDir()/local()/name();

            if (file(parentObjectPath))
            {
                return parentObjectPath;
            }
        }

        if (!dir(path))
        {
            word newInstancePath = time().findInstancePath(instant(instance()));

            if (newInstancePath.size())
            {
                fileName fName
                (
                    rootPath()/caseName()
                   /newInstancePath/db_.dbDir()/local()/name()
                );

                if (file(fName))
                {
                    return fName;
                }
            }
        }
    }

    return fileName::null;
}


Foam::Istream* Foam::IOobject::objectStream()
{
    fileName fName = filePath();

    if (fName != fileName::null)
    {
        IFstream* isPtr = new IFstream(fName);

        if (isPtr->good())
        {
            return isPtr;
        }
        else
        {
            delete isPtr;
            return NULL;
        }
    }
    else
    {
        return NULL;
    }
}


bool Foam::IOobject::headerOk()
{
    bool ok = true;

    Istream* isPtr = objectStream();

    // If the stream has failed return
    if (!isPtr)
    {
        if (objectRegistry::debug)
        {
            Info
                << "IOobject::headerOk() : "
                << "file " << objectPath() << " could not be opened"
                << endl;
        }

        ok = false;
    }
    else
    {
        // Try reading header
        if (!readHeader(*isPtr))
        {
            if (objectRegistry::debug)
            {
                IOWarningIn("IOobject::headerOk()", (*isPtr))
                    << "failed to read header of file " << objectPath()
                    << endl;
            }

            ok = false;
        }
    }

    delete isPtr;

    return ok;
}


void Foam::IOobject::setBad(const string& s)
{
    if (objState_ != GOOD)
    {
        FatalErrorIn("IOobject::setBad(const string&)")
            << "recurrent failure for object " << s
            << exit(FatalError);
    }

    if (error::level)
    {
        Info<< "IOobject::setBad(const string&) : "
            << "broken object " << s << info() << endl;
    }

    objState_ = BAD;
}


void Foam::IOobject::operator=(const IOobject& io)
{
    name_ = io.name_;
    headerClassName_ = io.headerClassName_;
    note_ = io.note_;
    instance_ = io.instance_;
    local_ = io.local_;
    rOpt_ = io.rOpt_;
    wOpt_ = io.wOpt_;
    objState_ = io.objState_;
}


// ************************************************************************* //
