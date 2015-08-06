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

#include "sphereRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(sphereRefinement, 0);
addToRunTimeSelectionTable(objectRefinement, sphereRefinement, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sphereRefinement::sphereRefinement()
:
    objectRefinement(),
    centre_(),
    radius_(-1.0)
{}

sphereRefinement::sphereRefinement
(
    const word& name,
    const scalar cellSize,
    const direction additionalRefLevels,
    const point& centre,
    const scalar radius
)
:
    objectRefinement(),
    centre_(centre),
    radius_(radius)
{
    setName(name);
    setCellSize(cellSize);
    setAdditionalRefinementLevels(additionalRefLevels);
}

sphereRefinement::sphereRefinement
(
    const word& name,
    const dictionary& dict
)
:
    objectRefinement(name, dict)
{
    this->operator=(dict);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool sphereRefinement::intersectsObject(const boundBox& bb) const
{
    const point& c = (bb.max() + bb.min()) / 2.0;
    
    if( magSqr(c - centre_) < sqr(radius_) )
        return true;
    
    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary sphereRefinement::dict(bool ignoreType) const
{
    dictionary dict;

    if( additionalRefinementLevels() == 0 && cellSize() >= 0.0 )
    {
        dict.add("cellSize", cellSize());
    }
    else
    {
        dict.add("additionalRefinementLevels", additionalRefinementLevels());
    }

    dict.add("type", type());

    dict.add("centre", centre_);
    dict.add("radius", radius_);

    return dict;
}

void sphereRefinement::write(Ostream& os) const
{
    os  << " type:   " << type()
        << " centre: " << centre_
        << " radius: " << radius_;
}

void sphereRefinement::writeDict(Ostream& os, bool subDict) const
{
    if( subDict )
    {
        os << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }
    
    if( additionalRefinementLevels() == 0 && cellSize() >= 0.0 )
    {
        os.writeKeyword("cellSize") << cellSize() << token::END_STATEMENT << nl;
    }
    else
    {
        os.writeKeyword("additionalRefinementLevels")
                << additionalRefinementLevels()
                << token::END_STATEMENT << nl;
    }

    // only write type for derived types
    if( type() != typeName_() )
    {
        os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
    }

    os.writeKeyword("centre") << centre_ << token::END_STATEMENT << nl;
    os.writeKeyword("radius") << radius_ << token::END_STATEMENT << nl;

    if( subDict )
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}

void sphereRefinement::operator=(const dictionary& d)
{
    // allow as embedded sub-dictionary "coordinateSystem"
    const dictionary& dict =
    (
        d.found(typeName_())
      ? d.subDict(typeName_())
      : d
    );

    // unspecified centre is (0 0 0)
    if( dict.found("centre") )
    {
        dict.lookup("centre") >> centre_;
    }
    else
    {
        FatalErrorIn
        (
            "void sphereRefinement::operator=(const dictionary& d)"
        ) << "Entry centre is not specified!" << exit(FatalError);
        centre_ = vector::zero;
    }

    // specify radius
    if( dict.found("radius") )
    {
        radius_ = readScalar(dict.lookup("radius"));
    }
    else
    {
        FatalErrorIn
        (
            "void sphereRefinement::operator=(const dictionary& d)"
        ) << "Entry radius is not specified!" << exit(FatalError);
        radius_ = -1.0;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& sphereRefinement::operator<<(Ostream& os) const
{
    os << "name " << name() << nl;
    os << "cell size " << cellSize() << nl;
    os << "additionalRefinementLevels " << additionalRefinementLevels() << endl;

    write(os);
    return os;
}
        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        
} // End namespace Foam

// ************************************************************************* //
