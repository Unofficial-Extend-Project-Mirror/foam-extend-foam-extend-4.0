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

#include "coordinateSystems.H"
#include "IOPtrList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<coordinateSystem>, 0);
}

const Foam::word Foam::coordinateSystems::dataType("coordinateSystem");

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateSystems::coordinateSystems()
{}


Foam::coordinateSystems::coordinateSystems
(
    const objectRegistry& registry,
    const word& name,
    const fileName& instance
)
{
    IOPtrList<coordinateSystem> newList
    (
        IOobject
        (
            name,
            instance,
            registry,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false               // don't register
        )
    );

    transfer(newList);
}


Foam::coordinateSystems::coordinateSystems
(
    const IOobject& io
)
{
    IOPtrList<coordinateSystem> newList(io);
    transfer(newList);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::label Foam::coordinateSystems::find(const word& keyword) const
{
    forAll(*this, i)
    {
        if (keyword == operator[](i).name())
        {
            return i;
        }
    }

    return -1;
}


bool Foam::coordinateSystems::found(const word& keyword) const
{
    return find(keyword) >= 0;
}


Foam::wordList Foam::coordinateSystems::toc() const
{
    wordList keywords(size());

    forAll(*this, i)
    {
        keywords[i] = operator[](i).name();
    }

    return keywords;
}


bool Foam::coordinateSystems::rewriteDict(dictionary& dict, bool noType) const
{
    if (dict.found(dataType) && !dict.isDict(dataType))
    {
        word name(dict.lookup(dataType));
        label i = find(name);

        if (i >= 0)
        {
            dict.remove(dataType);
            dict.add(dataType, operator[](i).dict(noType));
            return true;
        }

        FatalErrorIn
        (
            "Foam::coordinateSystems::rewriteDict(dictionary&, bool) const"
        )   << "could not rewrite " << dataType << " " << name << nl
            << "available coordinate systems: " << toc() << nl << nl
            << "context: " << nl
            << dict << nl
            << exit(FatalError);
    }

    return false;
}


bool Foam::coordinateSystems::writeData(Ostream& os, bool subDict) const
{
    // Write size of list
    os << nl << size();

    // Write beginning of contents
    os << nl << token::BEGIN_LIST;

    // Write list contents
    forAll(*this, i)
    {
        os << nl;
        operator[](i).writeDict(os, subDict);
    }

    // Write end of contents
    os << token::END_LIST << nl;

    // Check state of IOstream
    return os.good();
}


// ************************************************************************* //
