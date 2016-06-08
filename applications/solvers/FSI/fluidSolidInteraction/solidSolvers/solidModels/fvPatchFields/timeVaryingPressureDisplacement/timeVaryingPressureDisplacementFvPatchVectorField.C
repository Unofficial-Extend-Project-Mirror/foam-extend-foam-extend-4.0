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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "timeVaryingPressureDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
// #include "constitutiveModel.H"
// #include "solidSolver.H"
// #include "pRveUnsTotalLagrangianSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeVaryingPressureDisplacementFvPatchVectorField::
timeVaryingPressureDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    timeSeries_()
{}


timeVaryingPressureDisplacementFvPatchVectorField::
timeVaryingPressureDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF, dict),
    timeSeries_(dict)
{
    traction() = vector::zero;
    pressure() = timeSeries_(this->db().time().timeOutputValue());

    Info << "Creating time varying pressure displacement boundary conditions"
        << endl;
}


timeVaryingPressureDisplacementFvPatchVectorField::
timeVaryingPressureDisplacementFvPatchVectorField
(
    const timeVaryingPressureDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(tdpvf, p, iF, mapper),
    timeSeries_(tdpvf.timeSeries_)
{}


timeVaryingPressureDisplacementFvPatchVectorField::
timeVaryingPressureDisplacementFvPatchVectorField
(
    const timeVaryingPressureDisplacementFvPatchVectorField& tdpvf
)
:
    tractionDisplacementFvPatchVectorField(tdpvf),
    timeSeries_(tdpvf.timeSeries_)
{}


timeVaryingPressureDisplacementFvPatchVectorField::
timeVaryingPressureDisplacementFvPatchVectorField
(
    const timeVaryingPressureDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(tdpvf, iF),
    timeSeries_(tdpvf.timeSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeVaryingPressureDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    tractionDisplacementFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void timeVaryingPressureDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    tractionDisplacementFvPatchVectorField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void timeVaryingPressureDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    pressure() = timeSeries_(this->db().time().timeOutputValue());

    tractionDisplacementFvPatchVectorField::updateCoeffs();
}

// Write
void timeVaryingPressureDisplacementFvPatchVectorField::write
(
    Ostream& os
) const
{
    tractionDisplacementFvPatchVectorField::write(os);
    timeSeries_.write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    timeVaryingPressureDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
