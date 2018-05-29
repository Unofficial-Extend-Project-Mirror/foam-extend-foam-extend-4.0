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

\*---------------------------------------------------------------------------*/

#include "noSlipMovingWallFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noSlipMovingWallFvPatchVectorField::noSlipMovingWallFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


noSlipMovingWallFvPatchVectorField::noSlipMovingWallFvPatchVectorField
(
    const noSlipMovingWallFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


noSlipMovingWallFvPatchVectorField::noSlipMovingWallFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


noSlipMovingWallFvPatchVectorField::noSlipMovingWallFvPatchVectorField
(
    const noSlipMovingWallFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf)
{}


noSlipMovingWallFvPatchVectorField::noSlipMovingWallFvPatchVectorField
(
    const noSlipMovingWallFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noSlipMovingWallFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = dimensionedInternalField().mesh();

    if (mesh.changing())
    {
        const fvPatch& p = patch();
        const polyPatch& pp = p.patch();
        const pointField& oldPoints = mesh.oldPoints();

        vectorField oldFc(pp.size());

        forAll(oldFc, i)
        {
            oldFc[i] = pp[i].centre(oldPoints);
        }

        // Get wall-parallel mesh motion velocity from geometry
        vectorField Up =
            (pp.faceCentres() - oldFc)/mesh.time().deltaT().value();

        const volVectorField& U =
            mesh.lookupObject<volVectorField>
            (
                dimensionedInternalField().name()
            );

        scalarField phip =
            p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U));

        vectorField n = p.nf();
        const scalarField& magSf = p.magSf();
        scalarField Un = phip/(magSf + VSMALL);

        // Adjust for surface-normal mesh motion flux
        vectorField::operator=(Up + n*(Un - (n & Up)));
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


tmp<vectorField>
noSlipMovingWallFvPatchVectorField::gradientInternalCoeffs() const
{
    const vectorField nPatch = this->patch().nf();

    // Calculate the diagonal part of the transformation tensor
    vectorField implicitCoeffs(nPatch.size(), vector::zero);
    contractLinear(implicitCoeffs, I - nPatch*nPatch);

    return -this->patch().deltaCoeffs()*implicitCoeffs;
}


tmp<vectorField>
noSlipMovingWallFvPatchVectorField::gradientBoundaryCoeffs() const
{
    // Get necessary field
    const vectorField& UPatch = *this;
    const vectorField nPatch = this->patch().nf();
    const tensorField T = I - nPatch*nPatch;
    const vectorField UPatchInternal = this->patchInternalField();

    // Calculate the explicit contribution in the loop for performance
    vectorField explicitCoeffs(nPatch.size(), vector::zero);

    forAll(explicitCoeffs, faceI)
    {
        const tensor& curT = T[faceI];
        tensor curDiagT(tensor::zero);
        expandLinear(curDiagT, contractLinear(curT));

        explicitCoeffs[faceI] =
            (curT & UPatch[faceI])
          - (
                (curT - curDiagT)
              & UPatchInternal[faceI]
            );
    }

    return this->patch().deltaCoeffs()*explicitCoeffs;
}


void noSlipMovingWallFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    noSlipMovingWallFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
