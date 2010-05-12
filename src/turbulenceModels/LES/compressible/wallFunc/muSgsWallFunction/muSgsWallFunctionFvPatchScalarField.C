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

#include "muSgsWallFunctionFvPatchScalarField.H"
#include "LESModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

muSgsWallFunctionFvPatchScalarField::muSgsWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


muSgsWallFunctionFvPatchScalarField::muSgsWallFunctionFvPatchScalarField
(
    const muSgsWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


muSgsWallFunctionFvPatchScalarField::muSgsWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    Istream& is
)
:
    fixedValueFvPatchScalarField(p, iF, is)
{}


muSgsWallFunctionFvPatchScalarField::muSgsWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


muSgsWallFunctionFvPatchScalarField::muSgsWallFunctionFvPatchScalarField
(
    const muSgsWallFunctionFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{}


muSgsWallFunctionFvPatchScalarField::muSgsWallFunctionFvPatchScalarField
(
    const muSgsWallFunctionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void muSgsWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    const LESModel& sgsModel
        = db().lookupObject<LESModel>("LESProperties");

    scalar kappa = readScalar(sgsModel.lookup("kappa"));

    scalar E = readScalar(sgsModel.subDict("wallFunctionCoeffs").lookup("E"));

    const scalarField& ry = patch().deltaCoeffs();

    const fvPatchVectorField& U =
        patch().lookupPatchField<volVectorField, vector>("U");

    scalarField magUp = mag(U.patchInternalField() - U);

    const scalarField& muw =
        patch().lookupPatchField<volScalarField, scalar>("mu");

    const scalarField& rhow =
        patch().lookupPatchField<volScalarField, scalar>("rho");

    scalarField& muSgsw = *this;

    scalarField magFaceGradU = mag(U.snGrad());

    forAll(muSgsw, facei)
    {
        scalar magUpara = magUp[facei];

        scalar utau = sqrt
        (
            (muSgsw[facei] + muw[facei])
            *magFaceGradU[facei]/rhow[facei]
        );

        if(utau > 0)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = kappa*magUpara/utau;
                scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                scalar f =
                    - utau/(ry[facei]*muw[facei]/rhow[facei])
                    + magUpara/utau
                    + 1/E*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    - 1.0/(ry[facei]*muw[facei]/rhow[facei])
                    - magUpara/sqr(utau)
                    - 1/E*kUu*fkUu/utau;

                scalar utauNew = utau - f/df;
                err = mag((utau - utauNew)/utau);
                utau = utauNew;

            } while (utau > VSMALL && err > 0.01 && ++iter < 10);

            muSgsw[facei] =
                max(rhow[facei]*sqr(utau)/magFaceGradU[facei] - muw[facei],0.0);
        }
        else
        {
            muSgsw[facei] = 0;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, muSgsWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
