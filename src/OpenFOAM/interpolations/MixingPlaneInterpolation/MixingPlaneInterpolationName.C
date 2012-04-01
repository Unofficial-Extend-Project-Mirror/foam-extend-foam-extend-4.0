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

Description
    Shared template name for MixingPlane interpolation

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "MixingPlaneInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::MixingPlaneInterpolationName, 0);

template<>
const char*
Foam::NamedEnum
<
    Foam::MixingPlaneInterpolationName::assembly,
    4
>::names[] =
{
    "master",
    "slave",
    "both",
    "userDefined"
};


const Foam::NamedEnum
<
    Foam::MixingPlaneInterpolationName::assembly,
    4
>
Foam::MixingPlaneInterpolationName::assemblyNames_;


template<>
const char*
Foam::NamedEnum
<
    Foam::MixingPlaneInterpolationName::orientation,
    13
>::names[] =
{
    "unknown",
    "dirX_spanY",
    "dirX_spanZ",
    "dirY_spanX",
    "dirY_spanZ",
    "dirZ_spanX",
    "dirZ_spanY",

    "dirR_spanTheta",
    "dirR_spanZ",
    "dirTheta_spanR",
    "dirTheta_spanZ",
    "dirZ_spanR",
    "dirZ_spanTheta"
};


const Foam::NamedEnum
<
    Foam::MixingPlaneInterpolationName::orientation,
    13
>
Foam::MixingPlaneInterpolationName::orientationNames_;


template<>
const char*
Foam::NamedEnum
<
    Foam::MixingPlaneInterpolationName::mixingType,
    5
>::names[] =
{
    "averageFromNeighbourPatch",
    "averageFromNeighbourCellCenter",
    "averageFromOwnPatch",
    "zeroGradient",
    "doNothing"
};


const Foam::NamedEnum
<
    Foam::MixingPlaneInterpolationName::mixingType,
    5
>
Foam::MixingPlaneInterpolationName::mixingTypeNames_;

// ************************************************************************* //
