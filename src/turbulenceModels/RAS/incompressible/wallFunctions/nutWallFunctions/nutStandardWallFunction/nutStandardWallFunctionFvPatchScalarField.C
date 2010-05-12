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

#include "nutStandardWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutStandardWallFunctionFvPatchScalarField::
nutStandardWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("undefined")
{}


nutStandardWallFunctionFvPatchScalarField::
nutStandardWallFunctionFvPatchScalarField
(
    const nutStandardWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_)
{}


nutStandardWallFunctionFvPatchScalarField::
nutStandardWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{}


nutStandardWallFunctionFvPatchScalarField::
nutStandardWallFunctionFvPatchScalarField
(
    const nutStandardWallFunctionFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_)
{}


nutStandardWallFunctionFvPatchScalarField::
nutStandardWallFunctionFvPatchScalarField
(
    const nutStandardWallFunctionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nutStandardWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    const RASModel& rasModel
        = db().lookupObject<RASModel>("RASProperties");

    scalar kappa = rasModel.kappa().value();
    scalar E = rasModel.E().value();
    scalar yPlusLam = rasModel.yPlusLam();

    const scalarField& ry = patch().deltaCoeffs();

    const fvPatchVectorField& U =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    scalarField magUp = mag(U.patchInternalField() - U);

    const scalarField& nuw =
        patch().lookupPatchField<volScalarField, scalar>("nu");
    scalarField& nutw = *this;

    scalarField magFaceGradU = mag(U.snGrad());

    forAll(nutw, facei)
    {
        scalar magUpara = magUp[facei];

        scalar kappaRe = kappa*magUpara/(nuw[facei]*ry[facei]);

        scalar yPlus = yPlusLam;
        scalar ryPlusLam = 1.0/yPlus;

        int iter = 0;
        scalar yPlusLast = 0.0;

        do
        {
            yPlusLast = yPlus;
            yPlus = (kappaRe + yPlus)/(1.0 + log(E*yPlus));

        } while(mag(ryPlusLam*(yPlus - yPlusLast)) > 0.01 && ++iter < 10 );

        if (yPlus > yPlusLam)
        {
            nutw[facei] = nuw[facei]*(yPlus*kappa/log(E*yPlus) - 1);
        }
        else
        {
            nutw[facei] = 0.0;
        }
    }
}


void nutStandardWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
    os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutStandardWallFunctionFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
