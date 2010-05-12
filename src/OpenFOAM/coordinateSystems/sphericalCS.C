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

#include "sphericalCS.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sphericalCS, 0);
    addToRunTimeSelectionTable(coordinateSystem, sphericalCS, origAxisDir);
    addToRunTimeSelectionTable(coordinateSystem, sphericalCS, origRotation);
    addToRunTimeSelectionTable(coordinateSystem, sphericalCS, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sphericalCS::sphericalCS()
:
    coordinateSystem()
{}

Foam::sphericalCS::sphericalCS
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& direction
)
:
    coordinateSystem(name, origin, axis, direction)
{}


Foam::sphericalCS::sphericalCS
(
    const word& name,
    const point& origin,
    const coordinateRotation& cr
)
:
    coordinateSystem(name, origin, cr)
{}


Foam::sphericalCS::sphericalCS
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::vector Foam::sphericalCS::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    scalar r = local.x();
    scalar theta = local.y()*mathematicalConstant::pi/180.0;
    scalar phi = local.z()*mathematicalConstant::pi/180.0;

    return coordinateSystem::localToGlobal
    (
        vector(r*cos(theta)*sin(phi), r*sin(theta)*sin(phi), r*cos(phi)),
        translate
    );
}

Foam::tmp<Foam::vectorField> Foam::sphericalCS::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    const scalarField r = local.component(vector::X);

    const scalarField theta =
        local.component(vector::Y)*mathematicalConstant::pi/180.0;

    const scalarField phi =
        local.component(vector::Z)*mathematicalConstant::pi/180.0;

    vectorField lc(local.size());
    lc.replace(vector::X, r*cos(theta)*sin(phi));
    lc.replace(vector::Y, r*sin(theta)*sin(phi));
    lc.replace(vector::Z, r*cos(phi));

    return coordinateSystem::localToGlobal(lc, translate);
}


Foam::vector Foam::sphericalCS::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    const vector lc = coordinateSystem::globalToLocal(global, translate);
    const scalar r = mag(lc);

    return vector
    (
        r,
        atan2(lc.y(), lc.x())*180.0/mathematicalConstant::pi,
        acos(lc.z()/(r + SMALL))*180.0/mathematicalConstant::pi
    );
}


Foam::tmp<Foam::vectorField> Foam::sphericalCS::globalToLocal
(
    const vectorField& global,
    bool translate
) const
{
    const vectorField lc = coordinateSystem::globalToLocal(global, translate);
    const scalarField r = mag(lc);

    tmp<vectorField> tresult(new vectorField(lc.size()));
    vectorField& result = tresult();

    result.replace
    (
        vector::X, r

    );

    result.replace
    (
        vector::Y,
        atan2
        (
            lc.component(vector::Y),
            lc.component(vector::X)
        )*180.0/mathematicalConstant::pi
    );

    result.replace
    (
        vector::Z,
        acos(lc.component(vector::Z)/(r + SMALL))*180.0/mathematicalConstant::pi
    );

    return tresult;
}


// ************************************************************************* //
