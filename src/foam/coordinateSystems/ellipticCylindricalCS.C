/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ellipticCylindricalCS.H"

#include "Switch.H"
#include "mathematicalConstants.H"
#include "boundBox.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ellipticCylindricalCS, 0);
    addToRunTimeSelectionTable
    (
        coordinateSystem,
        ellipticCylindricalCS,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ellipticCylindricalCS::ellipticCylindricalCS(const bool inDegrees)
:
    coordinateSystem(),
    a_(0),
    inDegrees_(inDegrees)
{}


Foam::ellipticCylindricalCS::ellipticCylindricalCS
(
    const coordinateSystem& cs,
    const bool inDegrees
)
:
    coordinateSystem(cs),
    a_(0),
    inDegrees_(inDegrees)
{}


Foam::ellipticCylindricalCS::ellipticCylindricalCS
(
    const word& name,
    const coordinateSystem& cs,
    const bool inDegrees
)
:
    coordinateSystem(name, cs),
    a_(0),
    inDegrees_(inDegrees)
{}


Foam::ellipticCylindricalCS::ellipticCylindricalCS
(
    const word& name,
    const point& origin,
    const coordinateRotation& cr,
    const scalar a,
    const bool inDegrees
)
:
    coordinateSystem(name, origin, cr),
    a_(a),
    inDegrees_(inDegrees)
{}


Foam::ellipticCylindricalCS::ellipticCylindricalCS
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& direction,
    const scalar a,
    const bool inDegrees
)
:
    coordinateSystem(name, origin, axis, direction),
    a_(a),
    inDegrees_(inDegrees)
{}


Foam::ellipticCylindricalCS::ellipticCylindricalCS
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict),
    a_(readScalar(dict.lookup("a"))),
    inDegrees_(dict.lookupOrDefault<Switch>("degrees", true))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::coordinateSystem::spanInfo
Foam::ellipticCylindricalCS::spanLimited() const
{
    spanInfo b(Pair<bool>(true, true));

    // Upper bound or r is unlimited
    b[0] = Pair<bool>(true, false);

    // z is unlimited
    b[2] = Pair<bool>(false, false);

    return b;
}


Foam::boundBox Foam::ellipticCylindricalCS::spanBounds() const
{
    return boundBox
    (
        vector
        (
            0,
            ( inDegrees_ ? -180.0 : -mathematicalConstant::pi ),
            -GREAT
        ),
        vector
        (
            GREAT,
            ( inDegrees_ ? 180.0 : mathematicalConstant::pi ),
            GREAT
        )
    );
}


bool Foam::ellipticCylindricalCS::inDegrees() const
{
    return inDegrees_;
}


Foam::Switch& Foam::ellipticCylindricalCS::inDegrees()
{
    return inDegrees_;
}


Foam::vector Foam::ellipticCylindricalCS::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    // Notation: u = local.x() v = local.y() z = local.z();
    scalar theta =
        local.y()*( inDegrees_ ? mathematicalConstant::pi/180.0 : 1.0 );

    return coordinateSystem::localToGlobal
    (
        vector
        (
            a_*cosh(local.x())*cos(theta),
            a_*sinh(local.x())*sin(theta),
            local.z()
        ),
        translate
    );
}

Foam::tmp<Foam::vectorField> Foam::ellipticCylindricalCS::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    scalarField theta =
        local.component(vector::Y)*
        ( inDegrees_ ? mathematicalConstant::pi/180.0 : 1.0 );

    vectorField lc(local.size());
    lc.replace
    (
        vector::X,
        a_*cosh(local.component(vector::X))*cos(theta)
    );

    lc.replace
    (
        vector::Y,
        a_*sinh(local.component(vector::X))*sin(theta)
    );

    lc.replace
    (
        vector::Z,
        local.component(vector::Z)
    );

    return coordinateSystem::localToGlobal(lc, translate);
}


Foam::vector Foam::ellipticCylindricalCS::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    notImplemented
    (
        "ellipticCylindricalCS::globalToLocal(const vector&, bool) const"
    );

    return vector::zero;
}

Foam::tmp<Foam::vectorField> Foam::ellipticCylindricalCS::globalToLocal
(
    const vectorField& global,
    bool translate
) const
{
    notImplemented
    (
        "ellipticCylindricalCS::globalToLocal(const vectorField&, bool) const"
    );

    return tmp<vectorField>(vectorField::null());
}


void Foam::ellipticCylindricalCS::write(Ostream& os) const
{
    coordinateSystem::write(os);
    os << " inDegrees: " << inDegrees() << endl;
}


void Foam::ellipticCylindricalCS::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << indent << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    coordinateSystem::writeDict(os, false);
    os.writeKeyword("inDegrees") << inDegrees_ << token::END_STATEMENT << nl;

    if (subDict)
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}


// ************************************************************************* //
