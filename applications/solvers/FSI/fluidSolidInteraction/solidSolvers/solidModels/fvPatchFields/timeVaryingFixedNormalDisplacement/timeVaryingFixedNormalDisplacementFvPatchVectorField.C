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

#include "timeVaryingFixedNormalDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "componentMixedPointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeVaryingFixedNormalDisplacementFvPatchVectorField::
timeVaryingFixedNormalDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedNormalDisplacementFvPatchVectorField(p, iF),
    timeSeries_()
{}


timeVaryingFixedNormalDisplacementFvPatchVectorField::
timeVaryingFixedNormalDisplacementFvPatchVectorField
(
    const timeVaryingFixedNormalDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedNormalDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    timeSeries_(ptf.timeSeries_)
{}


timeVaryingFixedNormalDisplacementFvPatchVectorField::
timeVaryingFixedNormalDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedNormalDisplacementFvPatchVectorField(p, iF, dict),
    timeSeries_(dict)
{}


timeVaryingFixedNormalDisplacementFvPatchVectorField::
timeVaryingFixedNormalDisplacementFvPatchVectorField
(
    const timeVaryingFixedNormalDisplacementFvPatchVectorField& tvfdpvf
)
:
    fixedNormalDisplacementFvPatchVectorField(tvfdpvf),
    timeSeries_(tvfdpvf.timeSeries_)
{}


timeVaryingFixedNormalDisplacementFvPatchVectorField::
timeVaryingFixedNormalDisplacementFvPatchVectorField
(
    const timeVaryingFixedNormalDisplacementFvPatchVectorField& tvfdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedNormalDisplacementFvPatchVectorField(tvfdpvf, iF),
    timeSeries_(tvfdpvf.timeSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeVaryingFixedNormalDisplacementFvPatchVectorField::updateCoeffs()
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

//     if (fieldName() == "DD")
//     {
//         const fvPatchField<vector>& D =
//             patch().lookupPatchField<volVectorField, vector>("D");
//         disp -= U;
//     }
//     else if (fieldName() != "D")
//     {
//         FatalError << "The displacement field should be U or DU"
//             << exit(FatalError);
//     }


//     Info << disp << endl;

    this->refValue() = disp;


    // Update componentMixed boundary condition on point displacement field

    word fieldName = this->dimensionedInternalField().name();

    pointVectorField& pointField =
        const_cast<pointVectorField&>
        (
            this->db().lookupObject<pointVectorField>("point" + fieldName)
        );

    if
    (
        pointField.boundaryField()[this->patch().index()].type()
     == componentMixedPointPatchVectorField::typeName
    )
    {
//         Info << "Setting field point" << fieldName << " at " <<
//             this->patch().name() << endl;

        componentMixedPointPatchVectorField& patchPointField =
            refCast<componentMixedPointPatchVectorField>
            (
                pointField.boundaryField()[this->patch().index()]
            );

        vectorField pointDisp
        (
            this->patch().patch().nPoints(),
            timeSeries_(this->db().time().timeOutputValue())
        );

//         Info << pointDisp << endl;

        patchPointField.refValue() = pointDisp;
    }


    fixedNormalDisplacementFvPatchVectorField::updateCoeffs();
}

void timeVaryingFixedNormalDisplacementFvPatchVectorField::
write(Ostream& os) const
{
    fixedNormalDisplacementFvPatchVectorField::write(os);
    timeSeries_.write(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    timeVaryingFixedNormalDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
