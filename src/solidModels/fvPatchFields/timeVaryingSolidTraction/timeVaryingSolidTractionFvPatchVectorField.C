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

#include "timeVaryingSolidTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeVaryingSolidTractionFvPatchVectorField::
timeVaryingSolidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    timeSeries_()
{}


timeVaryingSolidTractionFvPatchVectorField::
timeVaryingSolidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
solidTractionFvPatchVectorField(p, iF),
timeSeries_(dict)
{
    fieldName() = dimensionedInternalField().name();
    traction() = vector::zero;
    pressure() = 0.0;

    nonLinear() =
        nonLinearGeometry::nonLinearNames_.read(dict.lookup("nonLinear"));

    //- the leastSquares has zero non-orthogonal correction
    //- on the boundary
    //- so the gradient scheme should be extendedLeastSquares
    if
    (
      Foam::word
      (
          dimensionedInternalField().mesh().schemesDict().gradScheme
          (
              "grad(" + fieldName() + ")"
          )
      ) != "extendedLeastSquares"
    )
    {
        Warning << "The gradScheme for " << fieldName()
            << " should be \"extendedLeastSquares 0\" for the boundary "
            << "non-orthogonal correction to be right" << endl;
    }
}


timeVaryingSolidTractionFvPatchVectorField::
timeVaryingSolidTractionFvPatchVectorField
(
    const timeVaryingSolidTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(stpvf, p, iF, mapper),
    timeSeries_(stpvf.timeSeries_)
{}


timeVaryingSolidTractionFvPatchVectorField::
timeVaryingSolidTractionFvPatchVectorField
(
    const timeVaryingSolidTractionFvPatchVectorField& stpvf
)
:
    solidTractionFvPatchVectorField(stpvf),
    timeSeries_(stpvf.timeSeries_)
{}


timeVaryingSolidTractionFvPatchVectorField::
timeVaryingSolidTractionFvPatchVectorField
(
    const timeVaryingSolidTractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(stpvf, iF),
    timeSeries_(stpvf.timeSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeVaryingSolidTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidTractionFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void timeVaryingSolidTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    solidTractionFvPatchVectorField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void timeVaryingSolidTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    traction() = timeSeries_(this->db().time().timeOutputValue());

    solidTractionFvPatchVectorField::updateCoeffs();
}

// Write
void timeVaryingSolidTractionFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);
    timeSeries_.write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    timeVaryingSolidTractionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
