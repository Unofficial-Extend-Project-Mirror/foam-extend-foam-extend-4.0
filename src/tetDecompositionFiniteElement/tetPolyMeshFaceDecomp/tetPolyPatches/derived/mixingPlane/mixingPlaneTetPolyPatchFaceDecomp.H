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
    mixingPlaneTetPolyPatchFaceDecomp

Description
    MixingPlane patch.

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    mixingPlaneTetPolyPatchFaceDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPlaneTetPolyPatchFaceDecomp_H
#define mixingPlaneTetPolyPatchFaceDecomp_H

#include "faceTetPolyPatchFaceDecomp.H"
#include "mixingPlanePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class mixingPlaneTetPolyPatchFaceDecomp Declaration
\*---------------------------------------------------------------------------*/

class mixingPlaneTetPolyPatchFaceDecomp
:
    public faceTetPolyPatchFaceDecomp
{

public:

    //- Runtime type information
    TypeName(mixingPlanePolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        mixingPlaneTetPolyPatchFaceDecomp
        (
            const polyPatch& patch,
            const tetPolyBoundaryMeshFaceDecomp& bm
        )
        :
            faceTetPolyPatchFaceDecomp(patch, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
