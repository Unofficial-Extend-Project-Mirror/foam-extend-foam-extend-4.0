/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved

Description
    Specilalisation of patchFlux member function for scalars needed to
    ensure conservation across GGI patches that are partially covered if the
    bridge overlap switch is on.

\*---------------------------------------------------------------------------*/

#include "ggiFvPatchField.H"
#include "surfaceFields.H"
#include "fvScalarMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<>
void Foam::ggiFvPatchField<scalar>::patchFlux
(
    GeometricField<scalar, fvsPatchField, surfaceMesh>& flux,
    const fvMatrix<scalar>& matrix
) const
{
    // Since we have adjusted the internal/boundary coefficients in the
    // manipulateMatrix member function below, we must not use
    // patchNeighbourField for reconstructing the flux. We only need to use
    // interpolated shadow field. VV, 6/Mar/2018.

    // Get patch ID
    const label patchI = this->patch().index();

    // Get internal field
    const scalarField& iField = this->internalField();

    // Get shadow face-cells and assemble shadow field
    const unallocLabelList& sfc = ggiPatch_.shadow().faceCells();

    scalarField sField(sfc.size());
    forAll (sField, i)
    {
        sField[i] = iField[sfc[i]];
    }

    // Interpolate shadow to this side. Note: must not bridge since internal
    // coeffs and boundary coeffs take it into account
    scalarField neighbourField(ggiPatch_.interpolate(sField));

    if (ggiPatch_.bridgeOverlap())
    {
        // Only set fully uncovered faces. Partially covered faces taken into
        // account by manipulating value and gradient matrix coefficients. Note:
        // mirror field is the same as patch internal field
        ggiPatch_.setUncoveredFaces
        (
            this->patchInternalField()(),
            neighbourField
        );
    }

    // Calculate the flux with correct neighbour field (fully uncovered faces
    // bridged, while partially uncovered faces taken into account by
    // manipulating value and gradient matrix coefficients in order to ensure
    // conservation for both convection and diffusion part across partially
    // overlapping faces). VV, 14/Mar/2018.
    flux.boundaryField()[patchI] =
        matrix.internalCoeffs()[patchI]*this->patchInternalField()
      - matrix.boundaryCoeffs()[patchI]*neighbourField;

    // Scale the flux on slave patch to ensure global conservation across this
    // partially overlapping GGI patch. The slight disbalance can happen since
    // we interpolate the matrix coeffs and field separately, and we should do
    // it together to ensure conservation. Current code design does not easily
    // allow this, so this is a meaningful temporary solution when we have
    // partially overlapping faces. VV, 14/Mar/2018.
    if (ggiPatch_.bridgeOverlap() && !ggiPatch_.master())
    {
        // This is slave patch, master already updated the fluxes and we can use
        // that info to scale the fluxes on this side

        // Get the total flux through master patch
        const scalar masterFlux =
            gSum(flux.boundaryField()[ggiPatch_.shadowIndex()]);

        // Get the total flux through slave patch
        const scalar slaveFlux = gSum(flux.boundaryField()[patchI]);

        // Calculate the scaling factor. Note: negative sign since master and
        // slave fluxes have opposite sign
        const scalar scalingFactor = -masterFlux/(slaveFlux + SMALL);

        // Scale the slave flux to ensure global patch conservation
        flux.boundaryField()[patchI] *= scalingFactor;

        if (debug)
        {
            Info<< "Scaling flux on patch: " << ggiPatch_.name()
                << " with " << scalingFactor
                << " to ensure conservation." << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
