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

#include "postfixedSubRegistry.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(postfixedSubRegistry, 0);
}


// * * * * * * * * * * * * * * * * Constructors *  * * * * * * * * * * * * * //

Foam::postfixedSubRegistry::postfixedSubRegistry
(
    const IOobject& io,
    const label nIoObjects
)
:
    objectRegistry(io, io.db().dbDir(), nIoObjects)
{
    writeOpt() = IOobject::NO_WRITE;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::postfixedSubRegistry::~postfixedSubRegistry()
{
    objectRegistry::clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::postfixedSubRegistry::mangleFileName(const fileName& fName) const
{
    return fileName(fName+name());
}


bool Foam::postfixedSubRegistry::checkIn(regIOobject& io) const
{
    if (postfixedSubRegistry::debug)
    {
        Pout<< "postfixedSubRegistry::checkIn(regIOobject&) : "
            << name() << " : checking in " << io.name()
            << endl;
    }

    word demangledName = io.name();
    demangledName = demangledName(io.name().size() - name().size());

    // Check into Mesh-Registry under full name
    const_cast<objectRegistry&>(parent()).insert(io.name(), &io);

    // Check into Mesh-Registry under short name
    return const_cast<postfixedSubRegistry&>(*this).insert(demangledName, &io);;
}


bool Foam::postfixedSubRegistry::checkOut(regIOobject& io) const
{
    word demangledName = io.name();
    demangledName = demangledName(io.name().size() - name().size());
    iterator iter = const_cast<postfixedSubRegistry&>(*this).find(demangledName);
    iterator iterPar = const_cast<objectRegistry&>(parent()).find(io.name());

    if (iter != end() && iterPar != end())
    {
        if (postfixedSubRegistry::debug)
        {
            Pout<< "postfixedSubRegistry::checkOut(regIOobject&) : "
                << name() << " : checking out " << io.name()
                << endl;
        }

        if (iter() != &io || iterPar() != &io)
        {
            if (postfixedSubRegistry::debug)
            {
                WarningIn("postfixedSubRegistry::checkOut(regIOobject&)")
                    << name() << " : attempt to checkOut copy of " << io.name()
                    << endl;
            }

            return false;
        }
        else
        {
            bool hasErased =
                const_cast<postfixedSubRegistry&>(*this).erase(iter)
                && const_cast<objectRegistry&>(parent()).erase(iterPar);

            if (io.ownedByRegistry())
            {
                delete iter();
            }

            return hasErased;
        }
    }
    else
    {
        if (postfixedSubRegistry::debug)
        {
            Pout<< "postfixedSubRegistry::checkOut(regIOobject&) : "
                << name() << " : could not find " << io.name()
                << " in registry " << name()
                << endl;
        }

        return false;
    }
}


// ************************************************************************* //
