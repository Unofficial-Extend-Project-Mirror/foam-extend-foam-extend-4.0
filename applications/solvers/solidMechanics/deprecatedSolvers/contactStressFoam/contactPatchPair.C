/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Description

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "contactPatchPair.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
    Foam::contactPatchPair::contactPatchPair
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


// ************************************************************************* //
