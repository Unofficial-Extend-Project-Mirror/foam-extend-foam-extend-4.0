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

#include "ellipticCylindricalCS.H"

#include "Switch.H"
#include "mathematicalConstants.H"
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


// ************************************************************************* //
