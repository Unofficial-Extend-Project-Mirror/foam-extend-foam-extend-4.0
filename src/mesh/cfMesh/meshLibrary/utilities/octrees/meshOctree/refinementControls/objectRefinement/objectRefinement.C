/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "IOstream.H"
#include "objectRefinement.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    
defineTypeNameAndDebug(objectRefinement, 0);
defineRunTimeSelectionTable(objectRefinement, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectRefinement::objectRefinement()
:
    name_(),
    cellSize_()
{}


objectRefinement::objectRefinement
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    cellSize_(readScalar(dict.lookup("cellSize")))
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

objectRefinement::~objectRefinement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
dictionary Foam::objectRefinement::dict(bool ignoreType) const
{
    dictionary dict;

    dict.add("name", name_);

    // only write type for derived types
    if (!ignoreType && type() != typeName_())
    {
        dict.add("type", type());
    }

    dict.add("origin", origin_);
    dict.add("e1", e1());
    dict.add("e3", e3());

    return dict;
}

void objectRefinement::write(Ostream& os) const
{
    os  << type()
        << " origin: " << origin()
        << " e1: " << e1() << " e3: " << e3();
}


void objectRefinement::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << indent << name_ << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    // only write type for derived types
    if (type() != typeName_())
    {
        os.writeKeyword("type")  << type()      << token::END_STATEMENT << nl;
    }

    os.writeKeyword("origin") << origin_  << token::END_STATEMENT << nl;
    os.writeKeyword("e1")     << e1()     << token::END_STATEMENT << nl;
    os.writeKeyword("e3")     << e3()     << token::END_STATEMENT << nl;

    if (subDict)
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}

void coordinateSystem::operator=(const dictionary& rhs)
{
    if (debug)
    {
        Pout<< "coordinateSystem::operator=(const dictionary&) : "
            << "assign from " << rhs << endl;
    }

    // allow as embedded sub-dictionary "coordinateSystem"
    const dictionary& dict =
    (
        rhs.found(typeName_())
      ? rhs.subDict(typeName_())
      : rhs
    );

    // unspecified origin is (0 0 0)
    if (dict.found("origin"))
    {
        dict.lookup("origin") >> origin_;
    }
    else
    {
        origin_ = vector::zero;
    }

    // specify via coordinateRotation
    if (dict.found("coordinateRotation"))
    {
        autoPtr<coordinateRotation> cr =
            coordinateRotation::New(dict.subDict("coordinateRotation"));

        R_  = cr();
    }
    else
    {
        // no sub-dictionary - specify via axes
        R_ = coordinateRotation(dict);
    }

    Rtr_ = R_.T();
}
*/

Ostream& operator<<(Ostream& os, const objectRefinement& obr)
{
    os << obr.name() << nl;
    obr.writeDict(os, true);
    return os;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
