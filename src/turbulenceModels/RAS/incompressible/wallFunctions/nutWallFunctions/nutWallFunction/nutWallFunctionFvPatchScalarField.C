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

#include "nutWallFunctionFvPatchScalarField.H"
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

nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("undefined")
{}


nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const nutWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_)
{}


nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{}


nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const nutWallFunctionFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_)
{}


nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const nutWallFunctionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nutWallFunctionFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    const RASModel& rasModel
        = db().lookupObject<RASModel>("RASProperties");

    scalar kappa = rasModel.kappa().value();
    scalar E = rasModel.E().value();

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

        scalar utau = sqrt((nutw[facei] + nuw[facei])*magFaceGradU[facei]);

        if (utau > VSMALL)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = min(kappa*magUpara/utau, 50);
                scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                scalar f =
                    - utau/(ry[facei]*nuw[facei])
                    + magUpara/utau
                    + 1/E*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    1.0/(ry[facei]*nuw[facei])
                  + magUpara/sqr(utau)
                  + 1/E*kUu*fkUu/utau;

                scalar utauNew = utau + f/df;
                err = mag((utau - utauNew)/utau);
                utau = utauNew;

            } while (utau > VSMALL && err > 0.01 && ++iter < 10);

            nutw[facei] =
                max(sqr(max(utau, 0))/magFaceGradU[facei] - nuw[facei], 0.0);
        }
        else
        {
            nutw[facei] = 0;
        }
    }
}


void nutWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
    os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, nutWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
