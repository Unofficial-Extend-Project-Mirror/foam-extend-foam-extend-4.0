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

#include "mutkWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> mutkWallFunctionFvPatchScalarField::calcMut() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalarField& y = turbModel.y()[patchI];
    const scalarField& rhow = turbModel.rho().boundaryField()[patchI];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const scalarField& muw = turbModel.mu().boundaryField()[patchI];

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tmutw(new scalarField(patch().size(), 0.0));
    scalarField& mutw = tmutw();

    // Get face cells
    const unallocLabelList& fc = patch().faceCells();

    forAll(mutw, faceI)
    {
        const label faceCellI = fc[faceI];

        const scalar yPlus =
            Cmu25*y[faceI]*sqrt(k[faceCellI])/(muw[faceI]/rhow[faceI]);

        if (yPlus > yPlusLam_)
        {
            mutw[faceI] = muw[faceI]*(yPlus*kappa_/log(E_*yPlus) - 1.0);
        }
    }

    return tmutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(p, iF)
{}


mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mutWallFunctionFvPatchScalarField(p, iF, dict)
{}


mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const mutkWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const mutkWallFunctionFvPatchScalarField& wfpsf
)
:
    mutWallFunctionFvPatchScalarField(wfpsf)
{}


mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const mutkWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> mutkWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalarField& y = turbModel.y()[patchI];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const scalarField kwc(k.boundaryField()[patchI].patchInternalField());
    const scalarField& muw = turbModel.mu().boundaryField()[patchI];
    const scalarField& rhow = turbModel.rho().boundaryField()[patchI];

    return pow025(Cmu_)*y*sqrt(kwc)/(muw/rhow);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, mutkWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
