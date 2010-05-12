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
    
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "contactPatchPair.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
contactPatchPair::contactPatchPair
(
    const fvMesh& m,
    const label master,
    const label slave,
    const scalar tolerance
)
:
    mesh_(m),
    masterPatchIndex_(master),
    slavePatchIndex_(slave),
    masterInterpolate_(m.boundaryMesh()[masterPatchIndex_]),
    slaveInterpolate_(m.boundaryMesh()[slavePatchIndex_]),
    patchToPatchInterpolate_
    (
        m.boundaryMesh()[masterPatchIndex_],
        m.boundaryMesh()[slavePatchIndex_],
        intersection::FULL_RAY,
        intersection::CONTACT_SPHERE
    ),
    tol_(tolerance),
    touchFraction_
    (
        mesh_.boundaryMesh()[slavePatchIndex_].size(),
        pTraits<scalar>::zero
    ),
    slaveDisplacement_
    (
        mesh_.boundaryMesh()[slavePatchIndex_].size(),
        pTraits<vector>::zero
    )
{}


} // End namespace Foam

// ************************************************************************* //
