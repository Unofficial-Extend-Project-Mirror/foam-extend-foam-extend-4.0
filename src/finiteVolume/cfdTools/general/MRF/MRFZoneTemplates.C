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
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector& origin = origin_.value();
    const vector rotVel = Omega();

    // Internal faces
    const vectorField& CfIn = Cf.internalField();
    const vectorField& SfIn = Sf.internalField();

    register label faceI, patchFaceI;

    forAll (internalFaces_, i)
    {
        faceI = internalFaces_[i];

        phi[faceI] -=
            rho[faceI]*(rotVel ^ (CfIn[faceI] - origin)) & SfIn[faceI];
    }

    // Included patches: reset the flux to exactly zero to avoid
    // round-off issues
    forAll (includedFaces_, patchI)
    {
        forAll (includedFaces_[patchI], i)
        {
            patchFaceI = includedFaces_[patchI][i];

            phi.boundaryField()[patchI][patchFaceI] = 0.0;
        }
    }

    // Excluded patches
    forAll (excludedFaces_, patchI)
    {
        forAll (excludedFaces_[patchI], i)
        {
            patchFaceI = excludedFaces_[patchI][i];

            phi.boundaryField()[patchI][patchFaceI] -=
                rho.boundaryField()[patchI][patchFaceI]
               *(rotVel ^ (Cf.boundaryField()[patchI][patchFaceI] - origin))
              & Sf.boundaryField()[patchI][patchFaceI];
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
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector& origin = origin_.value();
    const vector rotVel = Omega();

    // Internal faces
    const vectorField& CfIn = Cf.internalField();
    const vectorField& SfIn = Sf.internalField();

    register label faceI, patchFaceI;

    forAll (internalFaces_, i)
    {
        faceI = internalFaces_[i];
        
        phi[faceI] += (rotVel ^ (CfIn[faceI] - origin)) & SfIn[faceI];
    }

    // Included patches
    forAll (includedFaces_, patchI)
    {
        forAll (includedFaces_[patchI], i)
        {
            patchFaceI = includedFaces_[patchI][i];

            phi.boundaryField()[patchI][patchFaceI] +=
                (rotVel ^ (Cf.boundaryField()[patchI][patchFaceI] - origin))
              & Sf.boundaryField()[patchI][patchFaceI];
        }
    }

    // Excluded patches
    forAll (excludedFaces_, patchI)
    {
        forAll (excludedFaces_[patchI], i)
        {
            patchFaceI = excludedFaces_[patchI][i];

            phi.boundaryField()[patchI][patchFaceI] +=
                (rotVel ^ (Cf.boundaryField()[patchI][patchFaceI] - origin))
              & Sf.boundaryField()[patchI][patchFaceI];
        }
    }
}


// ************************************************************************* //
