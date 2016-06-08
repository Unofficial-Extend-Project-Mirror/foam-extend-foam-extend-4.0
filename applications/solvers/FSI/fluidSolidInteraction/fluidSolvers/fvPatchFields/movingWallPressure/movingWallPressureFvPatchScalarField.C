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

#include "movingWallPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingWallPressureFvPatchScalarField::
movingWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


movingWallPressureFvPatchScalarField::
movingWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF)
{
    fvPatchField<scalar>::operator=(patchInternalField());
}


movingWallPressureFvPatchScalarField::
movingWallPressureFvPatchScalarField
(
    const movingWallPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


movingWallPressureFvPatchScalarField::
movingWallPressureFvPatchScalarField
(
    const movingWallPressureFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf)
{}


movingWallPressureFvPatchScalarField::
movingWallPressureFvPatchScalarField
(
    const movingWallPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void movingWallPressureFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    word fieldName = dimensionedInternalField().name();

    const volScalarField& p =
        mesh.lookupObject<volScalarField>(fieldName);

    const fvPatchField<vector>& gradP =
        patch().lookupPatchField<volVectorField, vector>("grad(p)");

    const fvPatchField<vector>& ddtU =
        patch().lookupPatchField<volVectorField, vector>("ddt(U)");

//     vectorField delta = this->patch().delta();

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    scalarField dPP = (k&gradP.patchInternalField());

    scalarField nGradPb(this->patch().size(), 0);

    if (p.dimensions() == dimPressure)
    {
        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dimensionedScalar rho
        (
            transportProperties.lookup("rho")
        );

        nGradPb = -rho.value()*(n&ddtU);
//         nGradPb = -rho.value()*(n&prevAcceleration_);
    }
    else
    {
        nGradPb = -(n&ddtU);
//         nGradPb = -(n&prevAcceleration_);
    }


//     scalarField gradPb = -(n&ddtU);
//     scalarField gradPp = (n&gradP.patchInternalField());

    Info << "ddtUn, max: " << max(n&ddtU)
        << ", avg: " << average(n&ddtU)
        << ", min: " << min(n&ddtU) << endl;

    Field<scalar>::operator=
    (
//         this->patchInternalField()
        this->patchInternalField() + dPP
      + nGradPb/this->patch().deltaCoeffs()
//         + 0.5*(gradient() + gradPp)/this->patch().deltaCoeffs()
    );

    Info << "p, max: " << max(*this)
        << ", avg: " << average(*this)
        << ", min: " << min(*this) << endl;

    fvPatchField<scalar>::evaluate();
}


// void movingWallPressureFvPatchScalarField::updateCoeffs()
// {
//     if (updated())
//     {
//         return;
//     }

//     const fvMesh& mesh = this->patch().boundaryMesh().mesh();
//     word fieldName = dimensionedInternalField().name();
//     const volScalarField& p =
//         mesh.lookupObject<volScalarField>(fieldName);

//     const fvPatchField<vector>& ddtU =
//         patch().lookupPatchField<volVectorField, vector>("ddt(U)");

//     vectorField n = this->patch().nf();

//     if (p.dimensions() == dimPressure)
//     {
//         IOdictionary transportProperties
//         (
//             IOobject
//             (
//                 "transportProperties",
//                 mesh.time().constant(),
//                 mesh,
//                 IOobject::MUST_READ,
//                 IOobject::NO_WRITE
//             )
//         );

//         dimensionedScalar rho
//         (
//             transportProperties.lookup("rho")
//         );

//         gradient() = -rho.value()*(n&ddtU);
//     }
//     else
//     {
//         gradient() = -(n&ddtU);
//     }

//     fixedGradientFvPatchScalarField::updateCoeffs();
// }


void movingWallPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    movingWallPressureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
