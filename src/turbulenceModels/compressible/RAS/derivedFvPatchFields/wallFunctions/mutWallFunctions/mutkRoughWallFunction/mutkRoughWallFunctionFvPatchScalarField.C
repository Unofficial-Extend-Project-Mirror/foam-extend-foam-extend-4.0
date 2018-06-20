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

#include "mutkRoughWallFunctionFvPatchScalarField.H"
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

scalar mutkRoughWallFunctionFvPatchScalarField::fnRough
(
    const scalar KsPlus,
    const scalar Cs
) const
{
    // Return fn based on non-dimensional roughness height

    if (KsPlus < 90.0)
    {
        return pow
        (
            (KsPlus - 2.25)/87.75 + Cs*KsPlus,
            sin(0.4258*(log(KsPlus) - 0.811))
        );
    }
    else
    {
        return (1.0 + Cs*KsPlus);
    }
}


tmp<scalarField> mutkRoughWallFunctionFvPatchScalarField::calcMut() const
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

    tmp<scalarField> tmutw(new scalarField(*this));
    scalarField& mutw = tmutw();

    const unallocLabelList& fc = patch().faceCells();

    forAll(mutw, faceI)
    {
        const label faceCellI = fc[faceI];

        const scalar uStar = Cmu25*sqrt(k[faceCellI]);
        const scalar yPlus = uStar*y[faceI]/(muw[faceI]/rhow[faceI]);
        const scalar KsPlus = uStar*Ks_[faceI]/(muw[faceI]/rhow[faceI]);

        scalar Edash = E_;
        scalar yPlusLamNew = yPlusLam_;
        if (KsPlus > 2.25)
        {
            Edash /= fnRough(KsPlus, Cs_[faceI]);
            yPlusLamNew = turbModel.yPlusLam(kappa_, Edash);
        }

        if (yPlus > yPlusLamNew)
        {
            mutw[faceI] =
                muw[faceI]*(yPlus*kappa_/log(max(Edash*yPlus, 1 + 1e-4)) - 1);
        }

        if (debug)
        {
            Info<< "yPlus = " << yPlus
                << ", KsPlus = " << KsPlus
                << ", Edash = " << Edash
                << ", yPlusLam = " << yPlusLam_
                << ", mutw = " << mutw[faceI]
                << endl;
        }
    }

    return tmutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutkWallFunctionFvPatchScalarField(p, iF),
    Ks_(p.size(), 0.0),
    Cs_(p.size(), 0.0)
{}


mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const mutkRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Ks_(ptf.Ks_, mapper),
    Cs_(ptf.Cs_, mapper)
{}


mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mutkWallFunctionFvPatchScalarField(p, iF, dict),
    Ks_("Ks", dict, p.size()),
    Cs_("Cs", dict, p.size())
{}


mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const mutkRoughWallFunctionFvPatchScalarField& rwfpsf
)
:
    mutkWallFunctionFvPatchScalarField(rwfpsf),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_)
{}


mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const mutkRoughWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutkWallFunctionFvPatchScalarField(rwfpsf, iF),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mutkRoughWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mutkWallFunctionFvPatchScalarField::autoMap(m);
    Ks_.autoMap(m);
    Cs_.autoMap(m);
}


void mutkRoughWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const mutkRoughWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const mutkRoughWallFunctionFvPatchScalarField>(ptf);

    Ks_.rmap(nrwfpsf.Ks_, addr);
    Cs_.rmap(nrwfpsf.Cs_, addr);
}


void mutkRoughWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    Cs_.writeEntry("Cs", os);
    Ks_.writeEntry("Ks", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, mutkRoughWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
