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

#include "mutSpalartAllmarasWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "RASModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField>
mutSpalartAllmarasWallFunctionFvPatchScalarField::calcMut() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchI];
    const scalarField magGradU = mag(Uw.snGrad());
    const scalarField& rhow = turbModel.rho().boundaryField()[patchI];
    const scalarField& muw = turbModel.mu().boundaryField()[patchI];

    return Foam::max
    (
        scalar(0),
        rhow*sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - muw
    );
}


tmp<scalarField> mutSpalartAllmarasWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalarField& y = turbModel.y()[patchI];

    const fvPatchVectorField& Uw =
        turbModel.U().boundaryField()[patchI];
    const scalarField magUp = mag(Uw.patchInternalField() - Uw);

    const fvPatchScalarField& rhow = turbModel.rho().boundaryField()[patchI];

    const fvPatchScalarField& muw = turbModel.mu().boundaryField()[patchI];
    const scalarField& mutw = *this;

    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau();

    forAll(mutw, faceI)
    {
        const scalar& magUpara = magUp[faceI];

        scalar ut =
            sqrt((mutw[faceI] + muw[faceI])*magGradU[faceI]/rhow[faceI]);

        if (ut > VSMALL)
        {
            label iter = 0;
            scalar err = GREAT;

            do
            {
                const scalar kUu = min(kappa_*magUpara/ut, 50);
                const scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                const scalar f =
                    - ut*y[faceI]/(muw[faceI]/rhow[faceI])
                    + magUpara/ut
                    + 1/E_*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                const scalar df =
                    y[faceI]/(muw[faceI]/rhow[faceI])
                  + magUpara/sqr(ut)
                  + 1/E_*kUu*fkUu/ut;

                const scalar uTauNew = ut + f/df;
                err = mag((ut - uTauNew)/ut);
                ut = uTauNew;

            } while (ut > VSMALL && err > 0.01 && ++iter < 10);

            uTau[faceI] = max(0.0, ut);
        }
    }

    return tuTau;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutSpalartAllmarasWallFunctionFvPatchScalarField::
mutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(p, iF)
{}


mutSpalartAllmarasWallFunctionFvPatchScalarField::
mutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mutWallFunctionFvPatchScalarField(p, iF, dict)
{}


mutSpalartAllmarasWallFunctionFvPatchScalarField::
mutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const mutSpalartAllmarasWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


mutSpalartAllmarasWallFunctionFvPatchScalarField::
mutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const mutSpalartAllmarasWallFunctionFvPatchScalarField& wfpsf
)
:
    mutWallFunctionFvPatchScalarField(wfpsf)
{}


mutSpalartAllmarasWallFunctionFvPatchScalarField::
mutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const mutSpalartAllmarasWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField>
mutSpalartAllmarasWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalarField& y = turbModel.y()[patchI];
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchI];
    const scalarField& rhow = turbModel.rho().boundaryField()[patchI];
    const scalarField& muw = turbModel.mu().boundaryField()[patchI];

    return y*calcUTau(mag(Uw.snGrad()))/(muw/rhow);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    mutSpalartAllmarasWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
