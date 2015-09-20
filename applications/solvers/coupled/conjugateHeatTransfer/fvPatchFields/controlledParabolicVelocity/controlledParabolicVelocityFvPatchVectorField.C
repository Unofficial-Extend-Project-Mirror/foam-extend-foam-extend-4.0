/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "controlledParabolicVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

controlledParabolicVelocityFvPatchVectorField::controlledParabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Umean_(0),
    n_(0, 1, 0),
    y_(1, 0, 0),
    target_(0),
    obsFieldName_("undef"),
    obsPatchName_("undef"),
    obsPatchID_(-1),
    gain_(0),
    curTimeIndex_(-1)
{}


controlledParabolicVelocityFvPatchVectorField::controlledParabolicVelocityFvPatchVectorField
(
    const controlledParabolicVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Umean_(ptf.Umean_),
    n_(ptf.n_),
    y_(ptf.y_),
    target_(ptf.target_),
    obsFieldName_(ptf.obsFieldName_),
    obsPatchName_(ptf.obsPatchName_),
    obsPatchID_(ptf.obsPatchID_),
    gain_(ptf.gain_),
    curTimeIndex_(-1)
{}


controlledParabolicVelocityFvPatchVectorField::controlledParabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    Umean_(readScalar(dict.lookup("Umean"))),
    n_(dict.lookup("n")),
    y_(dict.lookup("y")),
    target_(readScalar(dict.lookup("target"))),
    obsFieldName_(dict.lookup("obsFieldName")),
    obsPatchName_(dict.lookup("obsPatchName")),
    obsPatchID_(patch().patch().boundaryMesh().findPatchID(obsPatchName_)),
    gain_(readScalar(dict.lookup("gain"))),
    curTimeIndex_(-1)
{
    if (obsPatchID_ < 0)
    {
        FatalErrorIn
        (
            "controlledParabolicVelocityFvPatchVectorField"
        )   << "patch " << dict.lookup("obsPatchName")
            << "not found"
            << nl << exit(FatalError);
    }

    n_ /= mag(n_);
    y_ /= mag(y_);

    fvPatchField<vector>::operator=
    (
        vectorField("value", dict, p.size())
    );
}


controlledParabolicVelocityFvPatchVectorField::controlledParabolicVelocityFvPatchVectorField
(
    const controlledParabolicVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    Umean_(fcvpvf.Umean_),
    n_(fcvpvf.n_),
    y_(fcvpvf.y_),
    target_(fcvpvf.target_),
    obsFieldName_(fcvpvf.obsFieldName_),
    obsPatchName_(fcvpvf.obsPatchName_),
    obsPatchID_(fcvpvf.obsPatchID_),
    gain_(fcvpvf.gain_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void controlledParabolicVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        const volScalarField& field =
            db().lookupObject<volScalarField>(obsFieldName_);

        scalar obsPatchAverage = gAverage(field.boundaryField()[obsPatchID_]);


        Info<< "Average of " << obsFieldName_
            << " on patch " << patch().name()
            << " = " << obsPatchAverage
            << " Difference to target = "
            << obsPatchAverage - target_
            << " Umean = " << Umean_
            << endl;

        Umean_ += gain_*(obsPatchAverage - target_);

        // Get range and orientation
        boundBox bb(patch().patch().localPoints(), true);

        vector ctr = 0.5*(bb.max() + bb.min());

        const vectorField& c = patch().Cf();

        // Calculate local 1-D coordinate for the parabolic profile
        scalarField coord =
            0.5 - ((c - ctr) & y_)/((bb.max() - bb.min()) & y_);

        operator==(n_*3/2*Umean_*(1.0 - sqr(coord)));

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void controlledParabolicVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("Umean")
        << Umean_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y")
        << y_ << token::END_STATEMENT << nl;
    os.writeKeyword("target")
        << target_ << token::END_STATEMENT << nl;
    os.writeKeyword("obsPatchName")
        << obsPatchName_ << token::END_STATEMENT << nl;
    os.writeKeyword("obsFieldName")
        << obsFieldName_ << token::END_STATEMENT << nl;
    os.writeKeyword("gain")
        << gain_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, controlledParabolicVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
