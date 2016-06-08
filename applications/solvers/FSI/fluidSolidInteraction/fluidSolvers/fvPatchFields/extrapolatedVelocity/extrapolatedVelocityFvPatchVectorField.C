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

\*---------------------------------------------------------------------------*/

#include "extrapolatedVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

extrapolatedVelocityFvPatchVectorField::
extrapolatedVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(p, iF)
{}


extrapolatedVelocityFvPatchVectorField::
extrapolatedVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchVectorField(p, iF)
{
    fvPatchField<vector>::operator=(patchInternalField());
}


extrapolatedVelocityFvPatchVectorField::
extrapolatedVelocityFvPatchVectorField
(
    const extrapolatedVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchVectorField(ptf, p, iF, mapper)
{}


extrapolatedVelocityFvPatchVectorField::
extrapolatedVelocityFvPatchVectorField
(
    const extrapolatedVelocityFvPatchVectorField& ptf
)
:
    zeroGradientFvPatchVectorField(ptf)
{}


extrapolatedVelocityFvPatchVectorField::
extrapolatedVelocityFvPatchVectorField
(
    const extrapolatedVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void extrapolatedVelocityFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>("grad(U)");

    vectorField delta = this->patch().delta();

    Field<vector>::operator=
    (
        this->patchInternalField() + (delta & gradU.patchInternalField())
    );

    fvPatchField<vector>::evaluate();
}


// void extrapolatedVelocityFvPatchVectorField::updateCoeffs()
// {
//     if (updated())
//     {
//         return;
//     }

//     const uniformDimensionedVectorField& g =
//         db().lookupObject<uniformDimensionedVectorField>("g");

//     const fvPatchField<vector>& rho =
//         patch().lookupPatchField<volVectorField, vector>(rhoName_);

//     // If the variable name is "p_rgh" or "pd" assume it is p - rho*g.h
//     // and set the gradient appropriately.
//     // Otherwise assume the variable is the static pressure.
//     if
//     (
//         dimensionedInternalField().name() == "p_rgh"
//      || dimensionedInternalField().name() == "pd"
//     )
//     {
//         gradient() = -rho.snGrad()*(g.value() & patch().Cf());
//     }
//     else
//     {
//         gradient() = rho*(g.value() & patch().nf());
//     }

//     zeroGradientFvPatchVectorField::updateCoeffs();
// }


void extrapolatedVelocityFvPatchVectorField::write(Ostream& os) const
{
    zeroGradientFvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    extrapolatedVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
