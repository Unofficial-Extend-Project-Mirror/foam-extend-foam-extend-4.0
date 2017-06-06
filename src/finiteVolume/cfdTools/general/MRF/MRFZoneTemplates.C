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

#include "MRFZone.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::MRFZone::relativeRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    // Get mesh velocity calculated from virtual mesh motion
    // HJ, 6/Jun/2017
    const surfaceScalarField& meshVel = meshVelocity();

    register label faceI, patchFaceI;

    // Internal faces
    scalarField& phiIn = phi.internalField();
    const scalarField& meshVelIn = meshVel.internalField();

    forAll (internalFaces_, i)
    {
        faceI = internalFaces_[i];

        phiIn[faceI] -= rho[faceI]*meshVelIn[faceI];
    }

    // Included patches: reset the flux to exactly zero to avoid
    // round-off issues
    forAll (includedFaces_, patchI)
    {
        forAll (includedFaces_[patchI], i)
        {
            patchFaceI = includedFaces_[patchI][i];

            phi.boundaryField()[patchI][patchFaceI] = 0;
        }
    }

    // Excluded patches
    forAll (excludedFaces_, patchI)
    {
        forAll (excludedFaces_[patchI], i)
        {
            patchFaceI = excludedFaces_[patchI][i];

            phi.boundaryField()[patchI][patchFaceI] -=
                rho.boundaryField()[patchI][patchFaceI]*
                meshVel.boundaryField()[patchI][patchFaceI];
        }
    }
}


template<class RhoFieldType>
void Foam::MRFZone::absoluteRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    // Get mesh velocity calculated from virtual mesh motion
    // HJ, 6/Jun/2017
    const surfaceScalarField& meshVel = meshVelocity();

    register label faceI, patchFaceI;

    // Internal faces
    scalarField& phiIn = phi.internalField();
    const scalarField& meshVelIn = meshVel.internalField();

    forAll (internalFaces_, i)
    {
        faceI = internalFaces_[i];

        phiIn[faceI] += meshVelIn[faceI];
    }

    // Included patches
    forAll (includedFaces_, patchI)
    {
        forAll (includedFaces_[patchI], i)
        {
            patchFaceI = includedFaces_[patchI][i];

            phi.boundaryField()[patchI][patchFaceI] +=
                meshVel.boundaryField()[patchI][patchFaceI];
        }
    }

    // Excluded patches
    forAll (excludedFaces_, patchI)
    {
        forAll (excludedFaces_[patchI], i)
        {
            patchFaceI = excludedFaces_[patchI][i];

            phi.boundaryField()[patchI][patchFaceI] +=
                meshVel.boundaryField()[patchI][patchFaceI];
        }
    }
}


// ************************************************************************* //
