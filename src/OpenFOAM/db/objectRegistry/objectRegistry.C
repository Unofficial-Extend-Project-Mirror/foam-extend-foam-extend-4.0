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

#include "objectRegistry.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(objectRegistry, 0);
}


// * * * * * * * * * * * * * * * * Constructors *  * * * * * * * * * * * * * //

Foam::objectRegistry::objectRegistry
(
    const Time& t,
    const label nIoObjects
)
:
    regIOobject
    (
        IOobject
        (
            string::validate<word>(t.caseName()),
            "",
            t,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        )
    ),
    HashTable<regIOobject*>(nIoObjects),
    time_(t),
    parent_(t),
    dbDir_(name())
{}


Foam::objectRegistry::objectRegistry
(
    const IOobject& io,
    const label nIoObjects
)
:
    regIOobject(io),
    HashTable<regIOobject*>(nIoObjects),
    time_(io.time()),
    parent_(io.db()),
    dbDir_(parent_.dbDir()/local()/name())
{
    writeOpt() = IOobject::AUTO_WRITE;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::objectRegistry::~objectRegistry()
{
    for (iterator iter = begin(); iter != end(); ++iter)
    {
        if (iter()->ownedByRegistry())
        {
            regIOobject* elemPtr = iter();
            erase(iter);
            delete elemPtr;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::objectRegistry::names() const
{
    wordList objectNames(size());

    label count=0;
    for (const_iterator iter = begin(); iter != end(); ++iter)
    {
        objectNames[count++] = iter()->name();
    }

    return objectNames;
}


Foam::wordList Foam::objectRegistry::names(const word& ClassName) const
{
    wordList objectNames(size());

    label count=0;
    for (const_iterator iter = begin(); iter != end(); ++iter)
    {
        if (iter()->type() == ClassName)
        {
            objectNames[count++] = iter()->name();
        }
    }

    objectNames.setSize(count);

    return objectNames;
}


const Foam::objectRegistry& Foam::objectRegistry::subRegistry
(
    const word& name
) const
{
    return lookupObject<objectRegistry>(name);
}


bool Foam::objectRegistry::checkIn(regIOobject& io) const
{
    if (objectRegistry::debug)
    {
        Pout<< "objectRegistry::checkIn(regIOobject&) : "
            << name() << " : checking in " << io.name()
            << endl;
    }

    return const_cast<objectRegistry&>(*this).insert(io.name(), &io);
}


bool Foam::objectRegistry::checkOut(regIOobject& io) const
{
    iterator iter = const_cast<objectRegistry&>(*this).find(io.name());

    if (iter != end())
    {
        if (objectRegistry::debug)
        {
            Pout<< "objectRegistry::checkOut(regIOobject&) : "
                << name() << " : checking out " << io.name()
                << endl;
        }

        if (iter() != &io)
        {
            if (objectRegistry::debug)
            {
                WarningIn("objectRegistry::checkOut(regIOobject&)")
                    << name() << " : attempt to checkOut copy of " << io.name()
                    << endl;
            }

            return false;
        }
        else
        {
            if (io.ownedByRegistry())
            {
                delete iter();
            }

            return const_cast<objectRegistry&>(*this).erase(iter);
        }
    }
    else
    {
        if (objectRegistry::debug)
        {
            Pout<< "objectRegistry::checkOut(regIOobject&) : "
                << name() << " : could not find " << io.name()
                << " in registry " << name()
                << endl;
        }

        return false;
    }
}


bool Foam::objectRegistry::modified() const
{
    bool anyModified = false;

    for (const_iterator iter = begin(); iter != end(); ++iter)
    {
        if (iter()->modified())
        {
            anyModified = true;
            break;
        }
    }

    return anyModified;
}


void Foam::objectRegistry::readModifiedObjects()
{
    for (iterator iter = begin(); iter != end(); ++iter)
    {
        if (objectRegistry::debug)
        {
            Pout<< "objectRegistry::readModifiedObjects() : "
                << name() << " : Considering reading object "
                << iter()->name()
                << endl;
        }

        iter()->readIfModified();
    }
}


bool Foam::objectRegistry::readIfModified()
{
    readModifiedObjects();
    return true;
}


bool Foam::objectRegistry::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    bool ok = true;

    for (const_iterator iter = begin(); iter != end(); ++iter)
    {
        if (objectRegistry::debug)
        {
            Pout<< "objectRegistry::write() : "
                << name() << " : Considering writing object "
                << iter()->name()
                << " with writeOpt " << iter()->writeOpt()
                << " to file " << iter()->objectPath()
                << endl;
        }

        if (iter()->writeOpt() != NO_WRITE)
        {
            ok = iter()->writeObject(fmt, ver, cmp) && ok;
        }
    }

    return ok;
}


// ************************************************************************* //
