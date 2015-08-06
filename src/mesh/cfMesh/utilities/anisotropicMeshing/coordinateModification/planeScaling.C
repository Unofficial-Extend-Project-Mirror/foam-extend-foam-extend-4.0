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

#include "planeScaling.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"
#include "plane.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(planeScaling, 0);
addToRunTimeSelectionTable
(
    coordinateModification,
    planeScaling,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

planeScaling::planeScaling()
:
    coordinateModification(),
    origin_(vector::zero),
    normal_(1, 1, 1),
    scalingDistance_(0.0),
    scalingFactor_(1.0)
{}

planeScaling::planeScaling
(
    const word& name,
    const point& origin,
    const vector& normal,
    const scalar scalingDistance,
    const scalar scalingFactor
)
:
    coordinateModification(),
    origin_(origin),
    normal_(normal/mag(normal)),
    scalingDistance_(scalingDistance),
    scalingFactor_(scalingFactor)
{
    if( scalingFactor_ < SMALL )
    {
        Warning << "Scaling factor for plane " << name << " is less than 0.0 "
                << endl;

        scalingFactor_= 1.0;
    }

    setName(name);
}

planeScaling::planeScaling
(
    const word& name,
    const dictionary& dict
)
:
    coordinateModification(name, dict)
{
    this->operator=(dict);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

point planeScaling::origin() const
{
    return origin_;
}

void planeScaling::translateAndModifyObject(const vector& disp)
{
    origin_ += disp;

    scalingDistance_ /= scalingFactor_;
}

vector planeScaling::displacement(const point& p) const
{
    const scalar dist = (p - origin_) & normal_;

    const vector translationVec =
        normal_ * scalingDistance_ * ((1.0/scalingFactor_) - 1.0);

    const scalar t = dist / scalingDistance_;

    const scalar tBnd = Foam::max(0.0, Foam::min(1.0, t));

    return tBnd * translationVec;
}

vector planeScaling::backwardDisplacement(const point& p) const
{
    const scalar dist = (p - origin_) & normal_;

    const vector translationVec =
        normal_ * scalingDistance_ * (scalingFactor_ - 1.0);

    const scalar t = dist / scalingDistance_;

    const scalar tBnd = Foam::max(0.0, Foam::min(1.0, t));

    return tBnd * translationVec;
}

bool planeScaling::combiningPossible() const
{
    return true;
}

void planeScaling::boundingPlanes(PtrList<plane>& pl) const
{
    if( Foam::mag(scalingFactor_ - 1.0) > VSMALL )
    {
        pl.setSize(2);

        pl.set(0, new plane(origin_, normal_));
        pl.set(1, new plane(origin_ + scalingDistance_ * normal_, normal_));
    }
    else
    {
        pl.clear();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary planeScaling::dict(bool ignoreType) const
{
    dictionary dict;

    dict.add("type", type());

    dict.add("origin", origin_);
    dict.add("normal", normal_);
    dict.add("scalingDistance", scalingDistance_);
    dict.add("scalingFactor", scalingFactor_);

    return dict;
}

void planeScaling::write(Ostream& os) const
{
    os  << " type:   " << type()
        << " origin: " << origin_
        << " normal: " << normal_
        << " scalingDistance: " << scalingDistance_
        << " scalingFactor: " << scalingFactor_;
}

void planeScaling::writeDict(Ostream& os, bool subDict) const
{
    if( subDict )
    {
        os << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    // only write type for derived types
    if( type() != typeName_() )
    {
        os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
    }

    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("normal") << normal_ << token::END_STATEMENT << nl;
    os.writeKeyword("scalingDistance") << scalingDistance_
                                           << token::END_STATEMENT << nl;
    os.writeKeyword("scalingFactor") << scalingFactor_
                                    << token::END_STATEMENT << nl;

    if( subDict )
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}

void planeScaling::operator=(const dictionary& d)
{
    // allow as embedded sub-dictionary "coordinateSystem"
    const dictionary& dict =
    (
        d.found(typeName_())
      ? d.subDict(typeName_())
      : d
    );

    // unspecified centre is (0 0 0)
    if( dict.found("origin") )
    {
        dict.lookup("origin") >> origin_;
    }
    else
    {
        FatalErrorIn
        (
            "void planeScaling::operator=(const dictionary& d)"
        ) << "Entry origin is not specified!" << exit(FatalError);

        origin_ = vector::zero;
    }

    // specify normal
    if( dict.found("normal") )
    {
        dict.lookup("normal") >> normal_;
    }
    else
    {
        FatalErrorIn
        (
            "void planeScaling::operator=(const dictionary& d)"
        ) << "Entry lengthX is not specified!" << exit(FatalError);

        normal_ = vector(1, 1, 1);
    }

    // specify translation distance
    if( dict.found("scalingDistance") )
    {
        scalingDistance_ = readScalar(dict.lookup("scalingDistance"));
    }
    else
    {
        FatalErrorIn
        (
            "void planeScaling::operator=(const dictionary& d)"
        ) << "Entry scalingDistance is not specified!" << exit(FatalError);

        scalingDistance_ = 0.0;
    }

    // specify scaling factor
    if( dict.found("scalingFactor") )
    {
        scalingFactor_ = readScalar(dict.lookup("scalingFactor"));
    }
    else
    {
        WarningIn
        (
            "void planeScaling::operator=(const dictionary& d)"
        ) << "Entry scalingFactor is not specified!" << endl;

        scalingFactor_ = 1.0;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& planeScaling::operator<<(Ostream& os) const
{
    os << "name " << name() << nl;
    write(os);
    return os;
}

Ostream& operator<<(Ostream& os, const planeScaling& pt)
{
    return pt.operator<<(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
