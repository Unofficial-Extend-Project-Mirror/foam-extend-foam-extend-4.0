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

#include "timeVaryingFixedDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeVaryingFixedDisplacementFvPatchVectorField::timeVaryingFixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    timeSeries_()
{}


timeVaryingFixedDisplacementFvPatchVectorField::timeVaryingFixedDisplacementFvPatchVectorField
(
    const timeVaryingFixedDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    timeSeries_(ptf.timeSeries_)
{}


timeVaryingFixedDisplacementFvPatchVectorField::timeVaryingFixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedDisplacementFvPatchVectorField(p, iF, dict),
    timeSeries_(dict)
{}


timeVaryingFixedDisplacementFvPatchVectorField::timeVaryingFixedDisplacementFvPatchVectorField
(
    const timeVaryingFixedDisplacementFvPatchVectorField& tvfdpvf
)
:
    fixedDisplacementFvPatchVectorField(tvfdpvf),
    timeSeries_(tvfdpvf.timeSeries_)
{}


timeVaryingFixedDisplacementFvPatchVectorField::timeVaryingFixedDisplacementFvPatchVectorField
(
    const timeVaryingFixedDisplacementFvPatchVectorField& tvfdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(tvfdpvf, iF),
    timeSeries_(tvfdpvf.timeSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeVaryingFixedDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    fvPatchField<vector>::operator==
    (
        timeSeries_(this->db().time().timeOutputValue())
    );
    fixedDisplacementFvPatchVectorField::updateCoeffs();
}

void timeVaryingFixedDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fixedDisplacementFvPatchVectorField::write(os);
    timeSeries_.write(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    timeVaryingFixedDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
