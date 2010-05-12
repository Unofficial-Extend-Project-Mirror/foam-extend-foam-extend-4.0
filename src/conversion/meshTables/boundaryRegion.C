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

#include "boundaryRegion.H"
#include "IOMap.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryRegion::boundaryRegion()
:
    Map<dictionary>()
{}


Foam::boundaryRegion::boundaryRegion
(
    const objectRegistry& registry,
    const word& name,
    const fileName& instance
)
:
    Map<dictionary>()
{
    readDict(registry, name, instance);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryRegion::~boundaryRegion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::boundaryRegion::append(const dictionary& dict)
{
    label maxId = -1;
    forAllConstIter(Map<dictionary>, *this, iter)
    {
        if (maxId < iter.key())
        {
            maxId = iter.key();
        }
    }

    insert(++maxId, dict);
    return maxId;
}


Foam::Map<Foam::word> Foam::boundaryRegion::names() const
{
    Map<word> lookup;

    forAllConstIter(Map<dictionary>, *this, iter)
    {
        word value = "boundaryRegion_" + Foam::name(iter.key());
        iter().readIfPresent("Label", value);

        lookup.insert(iter.key(), value);
    }

    return lookup;
}


Foam::Map<Foam::word> Foam::boundaryRegion::boundaryTypes() const
{
    Map<word> lookup;

    forAllConstIter(Map<dictionary>, *this, iter)
    {
        word value = "patch";
        iter().readIfPresent("BoundaryType", value);
        lookup.insert(iter.key(), value);
    }

    return lookup;
}


Foam::label Foam::boundaryRegion::findIndex(const word& name) const
{
    forAllConstIter(Map<dictionary>, *this, iter)
    {
        word theName;
        if (iter().readIfPresent("Label", theName))
        {
            if (theName == name)
            {
                return iter.key();
            }
        }
    }

    return -1;
}


Foam::word Foam::boundaryRegion::boundaryType(const word& name) const
{
    word bndType = "patch";

    label id = this->findIndex(name);
    if (id >= 0)
    {
        const dictionary& dict = operator[](id);
        dict.readIfPresent("BoundaryType", bndType);
    }

    return bndType;
}


void Foam::boundaryRegion::readDict
(
    const objectRegistry& registry,
    const word& name,
    const fileName& instance
)
{
    clear();

    // read constant/dictName
    IOMap<dictionary> ioObj
    (
        IOobject
        (
            name,
            instance,
            registry,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    if (ioObj.headerOk())
    {
        *this = ioObj;
    }
    else
    {
        Info<< "no constant/boundaryRegion information available" << endl;
    }
}


void Foam::boundaryRegion::writeDict
(
    const objectRegistry& registry,
    const word& name,
    const fileName& instance
) const
{
    // write constant/dictName
    IOMap<dictionary> ioObj
    (
        IOobject
        (
            name,
            instance,
            registry,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    ioObj.note() = "persistent data for thirdParty mesh <-> OpenFOAM translation";

    Info<< "Writing " << ioObj.name() << " to " << ioObj.objectPath() << endl;

    OFstream os(ioObj.objectPath());
    ioObj.writeHeader(os);
    os << *this;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::boundaryRegion::operator=(const boundaryRegion& rhs)
{
    Map<dictionary>::operator=(rhs);
}


void Foam::boundaryRegion::operator=(const Map<dictionary>& rhs)
{
    Map<dictionary>::operator=(rhs);
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

void Foam::boundaryRegion::rename(const dictionary& dict)
{
    if (!dict.size())
    {
        return;
    }

    Map<word> mapping;

    forAllConstIter(dictionary, dict, iter)
    {
        word oldName(iter().stream());

        label id = this->findIndex(oldName);
        if (id >= 0)
        {
            mapping.insert(id, iter().keyword());
        }
    }

    if (mapping.size())
    {
        forAllConstIter(Map<word>, mapping, iter)
        {
            label id = iter.key();
            word oldName(operator[](id).lookup("Label"));
            word newName(iter());

            dictionary newDict(operator[](id));
            newDict.remove("Label");
            newDict.add("Label", newName);

            this->set(id, newDict);
            Info<< "rename patch: " << newName << " <- " << oldName << endl;
        }
    }
}


// ************************************************************************* //
