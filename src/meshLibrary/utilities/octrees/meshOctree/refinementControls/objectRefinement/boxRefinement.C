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

#include "boxRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(boxRefinement, 0);
addToRunTimeSelectionTable(objectRefinement, boxRefinement, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

boxRefinement::boxRefinement()
:
    objectRefinement(),
    centre_(),
    lengthX_(-1.0),
    lengthY_(-1.0),
    lengthZ_(-1.0)
    //angleX_(0.0),
    //angleY_(0.0),
    //angleZ_(0.0)
{}

boxRefinement::boxRefinement
(
    const word& name,
    const scalar cellSize,
    const point& centre,
    const scalar lengthX,
    const scalar lengthY,
    const scalar lengthZ
    //const scalar angleX,
    //const scalar angleY,
    //const scalar angleZ
)
:
    objectRefinement(),
    centre_(centre),
    lengthX_(lengthX),
    lengthY_(lengthY),
    lengthZ_(lengthZ)
    //angleX_(angleX),
    //angleY_(angleY),
    //angleZ_(angleZ)
{
    setName(name);
    setCellSize(cellSize);
}

boxRefinement::boxRefinement
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

bool boxRefinement::intersectsObject(const boundBox& bb) const
{
    vector v(0.5*lengthX_, 0.5*lengthY_, 0.5*lengthZ_);
    boundBox box(centre_ - v, centre_ + v);
    
    if( box.overlaps(bb) )
        return true;
    
    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary boxRefinement::dict(bool ignoreType) const
{
    dictionary dict;

    dict.add("cellSize", cellSize());
    dict.add("type", type());

    dict.add("centre", centre_);
    dict.add("lengthX", lengthX_);
    dict.add("lengthY", lengthY_);
    dict.add("lengthZ", lengthZ_);
    
    //dict.add("angleX", angleX_);
    //dict.add("angleY", angleY_);
    //dict.add("angleZ", angleZ_);

    return dict;
}

void boxRefinement::write(Ostream& os) const
{
    os  << " type:   " << type()
        << " centre: " << centre_
        << " lengthX: " << lengthX_
        << " lengthY: " << lengthY_
        << " lengthZ: " << lengthZ_;
        //<< " angleX:  " << angleX_
        //<< " angleY:  " << angleY_
        //<< " angleZ:  " << angleZ_;
}

void boxRefinement::writeDict(Ostream& os, bool subDict) const
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

    os.writeKeyword("centre") << centre_ << token::END_STATEMENT << nl;
    os.writeKeyword("lengthX") << lengthX_ << token::END_STATEMENT << nl; 
    os.writeKeyword("lengthY") << lengthY_ << token::END_STATEMENT << nl;
    os.writeKeyword("lengthZ") << lengthZ_ << token::END_STATEMENT << nl;
    //os.writeKeyword("angleX") << angleX_ << token::END_STATEMENT << nl;
    //os.writeKeyword("angleY") << angleY_ << token::END_STATEMENT << nl;
    //os.writeKeyword("angleZ") << angleZ_ << token::END_STATEMENT << nl;

    if( subDict )
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}

void boxRefinement::operator=(const dictionary& d)
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
            "void boxRefinement::operator=(const dictionary& d)"
        ) << "Entry centre is not sopecified!" << exit(FatalError);
        centre_ = vector::zero;
    }

    // specify lengthX
    if( dict.found("lengthX") )
    {
        lengthX_ = readScalar(dict.lookup("lengthX"));
    }
    else
    {
        FatalErrorIn
        (
            "void boxRefinement::operator=(const dictionary& d)"
        ) << "Entry lengthX is not sopecified!" << exit(FatalError);
        lengthX_ = -1.0;
    }
    
    // specify lengthY
    if( dict.found("lengthY") )
    {
        lengthY_ = readScalar(dict.lookup("lengthY"));
    }
    else
    {
        FatalErrorIn
        (
            "void boxRefinement::operator=(const dictionary& d)"
        ) << "Entry lengthY is not sopecified!" << exit(FatalError);
        lengthY_ = -1.0;
    }
    
    // specify lengthZ
    if( dict.found("lengthZ") )
    {
        lengthZ_ = readScalar(dict.lookup("lengthZ"));
    }
    else
    {
        FatalErrorIn
        (
            "void boxRefinement::operator=(const dictionary& d)"
        ) << "Entry lengthZ is not sopecified!" << exit(FatalError);
        lengthZ_ = -1.0;
    }
    
    // specify angleX
/*    if( dict.found("angleX") )
    {
        angleX_ = readScalar(dict.lookup("angleX"));
    }
    else
    {
        angleX_ = -1.0;
    }
    
    // specify angleY
    if( dict.found("angleY") )
    {
        angleY_ = readScalar(dict.lookup("angleY"));
    }
    else
    {
        angleY_ = -1.0;
    }
    
    // specify angleX
    if( dict.found("angleZ") )
    {
        angleZ_ = readScalar(dict.lookup("angleZ"));
    }
    else
    {
        angleZ_ = -1.0;
    }
*/
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& boxRefinement::operator<<(Ostream& os) const
{
    os << "name " << name() << nl;
    os << "cell size " << cellSize() << nl;
    write(os);
    return os;
}
        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        
} // End namespace Foam

// ************************************************************************* //
