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

#include "coneRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(coneRefinement, 0);
addToRunTimeSelectionTable(objectRefinement, coneRefinement, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coneRefinement::coneRefinement()
:
    objectRefinement(),
    p0_(),
    r0_(-1.0),
    p1_(),
    r1_(-1.0)
{}

coneRefinement::coneRefinement
(
    const word& name,
    const scalar cellSize,
    const point& p0,
    const scalar radius0,
    const point& p1,
    const scalar radius1
)
:
    objectRefinement(),
    p0_(p0),
    r0_(radius0),
    p1_(p1),
    r1_(radius1)
{
    setName(name);
    setCellSize(cellSize);
}

coneRefinement::coneRefinement
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

bool coneRefinement::intersectsObject(const boundBox& bb) const
{
    //- check if the centre is inside the cone
    const point c = (bb.max() + bb.min()) / 2.0;
    
    const vector v = p1_ - p0_;
    const scalar d = magSqr(v);
    
    if( d < VSMALL )
        return false;
    
    const scalar t = ((c - p0_) & v) / d;
    if( (t > 1.0) || (t < 0.0) )
        return false;
    
    const scalar r = r0_ + (r1_ - r0_) * t;
    
    if( mag(p0_ + t * v - c) < r )
        return true;
    
    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary coneRefinement::dict(bool ignoreType) const
{
    dictionary dict;

    dict.add("cellSize", cellSize());
    dict.add("type", type());

    dict.add("p0", p0_);
    dict.add("radius0", r0_);
    dict.add("p1", p1_);
    dict.add("radius1", r1_);

    return dict;
}

void coneRefinement::write(Ostream& os) const
{
    os  << " type:   " << type()
        << " p0: " << p0_
        << " radius0: " << r0_
        << " p1: " << p1_
        << " radius1: " << r1_;
}

void coneRefinement::writeDict(Ostream& os, bool subDict) const
{
    if( subDict )
    {
        os << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }
    
    os.writeKeyword("cellSize") << cellSize() << token::END_STATEMENT << nl;

    // only write type for derived types
    if( type() != typeName_() )
    {
        os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
    }

    os.writeKeyword("p0") << p0_ << token::END_STATEMENT << nl;
    os.writeKeyword("radius0") << r0_ << token::END_STATEMENT << nl;
    os.writeKeyword("p1") << p1_ << token::END_STATEMENT << nl;
    os.writeKeyword("radius1") << r1_ << token::END_STATEMENT << nl;
    
    if( subDict )
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}

void coneRefinement::operator=(const dictionary& d)
{
    // allow as embedded sub-dictionary "coordinateSystem"
    const dictionary& dict =
    (
        d.found(typeName_())
      ? d.subDict(typeName_())
      : d
    );

    // unspecified centre is (0 0 0)
    if( dict.found("p0") )
    {
        dict.lookup("p0") >> p0_;
    }
    else
    {
        FatalErrorIn
        (
            "void coneRefinement::operator=(const dictionary& d)"
        ) << "Entry p0 is not specified!" << exit(FatalError);
        p0_ = vector::zero;
    }

    // specify radius
    if( dict.found("radius0") )
    {
        r0_ = readScalar(dict.lookup("radius0"));
    }
    else
    {
        FatalErrorIn
        (
            "void coneRefinement::operator=(const dictionary& d)"
        ) << "Entry radius0 is not specified!" << exit(FatalError);
        r0_ = -1.0;
    }
    
    // unspecified centre is (0 0 0)
    if( dict.found("p1") )
    {
        dict.lookup("p1") >> p1_;
    }
    else
    {
        FatalErrorIn
        (
            "void coneRefinement::operator=(const dictionary& d)"
        ) << "Entry p1 is not specified!" << exit(FatalError);
        p1_ = vector::zero;
    }

    // specify radius
    if( dict.found("radius1") )
    {
        r1_ = readScalar(dict.lookup("radius1"));
    }
    else
    {
        FatalErrorIn
        (
            "void coneRefinement::operator=(const dictionary& d)"
        ) << "Entry radius1 is not specified!" << exit(FatalError);
        r1_ = -1.0;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& coneRefinement::operator<<(Ostream& os) const
{
    os << "name " << name() << nl;
    os << "cell size " << cellSize() << nl;
    write(os);
    return os;
}
        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        
} // End namespace Foam

// ************************************************************************* //
