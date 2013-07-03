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

Class
    mixingPlaneInterpolation

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

Description
    Mixing plane interpolation: typedefs for polyPatch and faceZone

\*---------------------------------------------------------------------------*/

#ifndef mixingPlaneInterpolation_H
#define mixingPlaneInterpolation_H

#include "MixingPlaneInterpolation.H"
#include "PrimitivePatch.H"
#include "face.H"
#include "SubList.H"
#include "pointField.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    typedef MixingPlaneInterpolation
    <
        PrimitivePatch<face, SubList, const pointField&>,
        PrimitivePatch<face, SubList, const pointField&>
    >   mixingPlaneInterpolation;

    typedef MixingPlaneInterpolation
    <
        PrimitivePatch<face, List, const pointField&>,
        PrimitivePatch<face, List, const pointField&>
    >   mixingPlaneZoneInterpolation;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
