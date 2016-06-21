/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "timeVaryingFixedDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeVaryingFixedDisplacementFvPatchVectorField::
timeVaryingFixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    timeSeries_()
{}


timeVaryingFixedDisplacementFvPatchVectorField::
timeVaryingFixedDisplacementFvPatchVectorField
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


timeVaryingFixedDisplacementFvPatchVectorField::
timeVaryingFixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedDisplacementFvPatchVectorField(p, iF, dict),
    timeSeries_(dict)
{}


timeVaryingFixedDisplacementFvPatchVectorField::
timeVaryingFixedDisplacementFvPatchVectorField
(
    const timeVaryingFixedDisplacementFvPatchVectorField& tvfdpvf
)
:
    fixedDisplacementFvPatchVectorField(tvfdpvf),
    timeSeries_(tvfdpvf.timeSeries_)
{}


timeVaryingFixedDisplacementFvPatchVectorField::
timeVaryingFixedDisplacementFvPatchVectorField
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

    vectorField disp
        (
            patch().size(),
            timeSeries_(this->db().time().timeOutputValue())
            );

    if (fieldName() == "DU")
    {
        const fvPatchField<vector>& U =
          patch().lookupPatchField<volVectorField, vector>("U");
        disp -= U;
      }
    else if (fieldName() != "U")
    {
        FatalError
            << "The displacement field should be U or DU"
            << exit(FatalError);
    }

    fvPatchField<vector>::operator==
    (
     disp
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
