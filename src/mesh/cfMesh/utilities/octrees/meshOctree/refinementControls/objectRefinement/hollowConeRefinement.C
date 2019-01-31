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

#include "hollowConeRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(hollowConeRefinement, 0);
addToRunTimeSelectionTable(objectRefinement, hollowConeRefinement, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hollowConeRefinement::hollowConeRefinement()
:
    objectRefinement(),
    p0_(),
    r0Outer_(-1.0),
    r0Inner_(0.0),
    p1_(),
    r1Outer_(-1.0),
    r1Inner_(0.0)
{}

hollowConeRefinement::hollowConeRefinement
(
    const word& name,
    const scalar cellSize,
    const direction additionalRefLevels,
    const point& p0,
    const scalar radius0Outer,
    const scalar radius0Inner,
    const point& p1,
    const scalar radius1Outer,
    const scalar radius1Inner
)
:
    objectRefinement(),
    p0_(p0),
    r0Outer_(radius0Outer),
    r0Inner_(radius0Inner),
    p1_(p1),
    r1Outer_(radius1Outer),
    r1Inner_(radius1Inner)
{
    r0Inner_ = min(r0Inner_, r0Outer_);
    r1Inner_ = min(r0Inner_, r0Outer_);
    setName(name);
    setCellSize(cellSize);
    setAdditionalRefinementLevels(additionalRefLevels);
}

hollowConeRefinement::hollowConeRefinement
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

bool hollowConeRefinement::intersectsObject(const boundBox& bb) const
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
    
    const scalar rOuter = r0Outer_ + (r1Outer_ - r0Outer_) * t;
    const scalar rInner = r0Inner_ + (r1Inner_ - r0Inner_) * t;
    
    if(( mag(p0_ + t * v - c) < rOuter ) && ( mag(p0_ + t * v - c) > rInner ))
        return true;
    
    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary hollowConeRefinement::dict(bool /*ignoreType*/) const
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

    dict.add("p0", p0_);
    dict.add("radius0_Outer", r0Outer_);
    dict.add("radius0_Inner", r0Inner_);
    dict.add("p1", p1_);
    dict.add("radius1_Outer", r1Outer_);
    dict.add("radius1_Inner", r1Inner_);

    return dict;
}

void hollowConeRefinement::write(Ostream& os) const
{
    os  << " type:   " << type()
        << " p0: " << p0_
        << " radius0_Outer: " << r0Outer_
        << " radius0_Inner: " << r0Inner_
        << " p1: " << p1_
        << " radius1_Outer: " << r1Outer_
        << " radius1_Inner: " << r1Inner_;
}

void hollowConeRefinement::writeDict(Ostream& os, bool subDict) const
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

    os.writeKeyword("p0") << p0_ << token::END_STATEMENT << nl;
    os.writeKeyword("radius0_Outer") << r0Outer_ << token::END_STATEMENT << nl;
    os.writeKeyword("radius0_Inner") << r0Inner_ << token::END_STATEMENT << nl;
    os.writeKeyword("p1") << p1_ << token::END_STATEMENT << nl;
    os.writeKeyword("radius1_Outer") << r1Outer_ << token::END_STATEMENT << nl;
    os.writeKeyword("radius1_Inner") << r1Inner_ << token::END_STATEMENT << nl;
    
    if( subDict )
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}

void hollowConeRefinement::operator=(const dictionary& d)
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
            "void hollowConeRefinement::operator=(const dictionary& d)"
        ) << "Entry p0 is not specified!" << exit(FatalError);
        p0_ = vector::zero;
    }

    // specify radius
    if( dict.found("radius0_Outer") )
    {
        r0Outer_ = readScalar(dict.lookup("radius0_Outer"));
    }
    else
    {
        FatalErrorIn
        (
            "void hollowConeRefinement::operator=(const dictionary& d)"
        ) << "Entry radius0_Outer is not specified!" << exit(FatalError);
        r0Outer_ = -1.0;
    }

    // specify radius
    if( dict.found("radius0_Inner") )
    {
        r0Inner_ = readScalar(dict.lookup("radius0_Inner"));
    }
    else
    {
        FatalErrorIn
        (
     	    "void hollowConeRefinement::operator=(const dictionary& d)"
        ) << "Entry radius0_Inner is not specified!" << exit(FatalError);
        r0Inner_ = -1.0;
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
            "void hollowConeRefinement::operator=(const dictionary& d)"
        ) << "Entry p1 is not specified!" << exit(FatalError);
        p1_ = vector::zero;
    }

    // specify radius
    if( dict.found("radius1_Outer") )
    {
        r1Outer_ = readScalar(dict.lookup("radius1_Outer"));
    }
    else
    {
        FatalErrorIn
        (
            "void hollowConeRefinement::operator=(const dictionary& d)"
        ) << "Entry radius1_Outer is not specified!" << exit(FatalError);
        r1Outer_ = -1.0;
    }

    // specify radius
    if( dict.found("radius1_Inner") )
    {
     	r1Inner_ = readScalar(dict.lookup("radius1_Inner"));
    }
    else
    {
     	FatalErrorIn
        (
            "void hollowConeRefinement::operator=(const dictionary& d)"
        ) << "Entry radius1_Inner is not specified!" << exit(FatalError);
        r1Inner_ = -1.0;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& hollowConeRefinement::operator<<(Ostream& os) const
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
