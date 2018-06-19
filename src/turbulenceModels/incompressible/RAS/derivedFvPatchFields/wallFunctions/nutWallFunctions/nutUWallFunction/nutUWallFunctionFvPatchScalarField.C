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

#include "nutUWallFunctionFvPatchScalarField.H"
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutUWallFunctionFvPatchScalarField::calcYPlus
(
    const scalarField& magUp
) const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalarField& y = turbModel.y()[patchI];
    const scalarField& nuw = turbModel.nu().boundaryField()[patchI];

    tmp<scalarField> tyPlus(new scalarField(patch().size(), 0.0));
    scalarField& yPlus = tyPlus();

    forAll(yPlus, faceI)
    {
        const scalar kappaRe = kappa_*magUp[faceI]*y[faceI]/nuw[faceI];

        scalar yp = yPlusLam_;
        const scalar ryPlusLam = 1.0/yp;

        label iter = 0;
        scalar yPlusLast = 0.0;

        do
        {
            yPlusLast = yp;
            yp = (kappaRe + yp)/(1.0 + log(E_*yp));

        } while (mag(ryPlusLam*(yp - yPlusLast)) > 0.01 && ++iter < 10);

        yPlus[faceI] = max(0.0, yp);
    }

    return tyPlus;
}


tmp<scalarField> nutUWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchI];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    const scalarField& nuw = turbModel.nu().boundaryField()[patchI];

    tmp<scalarField> tyPlus = calcYPlus(magUp);
    scalarField& yPlus = tyPlus();

    tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
    scalarField& nutw = tnutw();

    forAll(yPlus, faceI)
    {
        if (yPlus[faceI] > yPlusLam_)
        {
            nutw[faceI] =
                nuw[faceI]*(yPlus[faceI]*kappa_/log(E_*yPlus[faceI]) - 1.0);
        }
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF)
{}


nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict)
{}


nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const nutUWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const nutUWallFunctionFvPatchScalarField& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf)
{}


nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const nutUWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutUWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchI = patch().index();
    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchI];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    return calcYPlus(magUp);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, nutUWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
