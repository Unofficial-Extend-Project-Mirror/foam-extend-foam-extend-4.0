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

#include "noSlipWallFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

noSlipWallFvPatchVectorField::noSlipWallFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


noSlipWallFvPatchVectorField::noSlipWallFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


noSlipWallFvPatchVectorField::noSlipWallFvPatchVectorField
(
    const noSlipWallFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


noSlipWallFvPatchVectorField::noSlipWallFvPatchVectorField
(
    const noSlipWallFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf)
{}


noSlipWallFvPatchVectorField::noSlipWallFvPatchVectorField
(
    const noSlipWallFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> noSlipWallFvPatchVectorField::gradientInternalCoeffs() const
{
    const vectorField nPatch = this->patch().nf();

    // Calculate the diagonal part of the transformation tensor
    vectorField implicitCoeffs(nPatch.size(), vector::zero);
    contractLinear(implicitCoeffs, I - nPatch*nPatch);

    return -this->patch().deltaCoeffs()*implicitCoeffs;
}


tmp<vectorField> noSlipWallFvPatchVectorField::gradientBoundaryCoeffs() const
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    noSlipWallFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
