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

#include "correctedSymmetryFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePatchFields(correctedSymmetry);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Evaluate the field on the patch
template<>
void correctedSymmetryFvPatchField<scalar>::evaluate(const Pstream::commsTypes)
{
    label secondOrder_ = false;

    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField nHat = this->patch().nf();

    vectorField delta = patch().delta();
    vectorField k = delta - nHat*(nHat&delta);

    word phiName = this->dimensionedInternalField().name();

    const fvPatchField<vector>& gradPhi =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + phiName + ")"
        );

    scalarField phiP = this->patchInternalField();
    phiP += (k & gradPhi.patchInternalField());

    if (secondOrder_)
    {
        scalarField nGradPhiP = (nHat & gradPhi.patchInternalField());

        Field<scalar>::operator=
        (
            phiP + 0.5*nGradPhiP/this->patch().deltaCoeffs()
        );
    }
    else
    {
        Field<scalar>::operator=(phiP);
    }

    transformFvPatchField<scalar>::evaluate();
}


// return gradient at boundary
template<>
tmp<Field<vector> > correctedSymmetryFvPatchField<vector>::snGrad() const
{
    label secondOrder_ = false;

    vectorField nHat = this->patch().nf();

    vectorField delta = patch().delta();
    vectorField k = delta - nHat*(nHat&delta);

    word phiName = this->dimensionedInternalField().name();

    const fvPatchField<tensor>& gradPhi =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + phiName + ")"
        );

    vectorField phiP = this->patchInternalField();
    phiP += (k & gradPhi.patchInternalField());

    if (secondOrder_)
    {
        vectorField nGradPhiP = (nHat & gradPhi.patchInternalField());

        return
          2*(
                transform(I - 2.0*sqr(nHat), phiP) - phiP
            )*(this->patch().deltaCoeffs()/2.0)
          - transform(sqr(nHat), nGradPhiP);
    }
    else
    {
        return
        (
            transform(I - 2.0*sqr(nHat), phiP)
          - phiP
        )*(this->patch().deltaCoeffs()/2.0);
    }
}


// Evaluate the field on the patch
template<>
void correctedSymmetryFvPatchField<vector>::evaluate(const Pstream::commsTypes)
{
    label secondOrder_ = false;

    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField nHat = this->patch().nf();

    vectorField delta = patch().delta();
    vectorField k = delta - nHat*(nHat&delta);

    word phiName = this->dimensionedInternalField().name();

    const fvPatchField<tensor>& gradPhi =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + phiName + ")"
        );

    vectorField phiP = this->patchInternalField();
    phiP += (k & gradPhi.patchInternalField());

    if (secondOrder_)
    {
        vectorField nGradPhiP = (nHat & gradPhi.patchInternalField());

        Field<vector>::operator=
        (
            transform
            (
                I - sqr(nHat),
                phiP + 0.5*nGradPhiP/this->patch().deltaCoeffs()
            )
        );
    }
    else
    {
        Field<vector>::operator=
        (
            (
                phiP
              + transform(I - 2.0*sqr(nHat), phiP)
            )/2.0
        );
    }

    transformFvPatchField<vector>::evaluate();
}



} // End namespace Foam

// ************************************************************************* //
