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

#include "SRFSurfaceNormalVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "SRFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SRFSurfaceNormalVelocityFvPatchVectorField::
SRFSurfaceNormalVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    relative_(0),
    refValue_(p.size(), 0)
{}


SRFSurfaceNormalVelocityFvPatchVectorField::
SRFSurfaceNormalVelocityFvPatchVectorField
(
    const SRFSurfaceNormalVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    relative_(ptf.relative_),
    refValue_(ptf.refValue_, mapper)
{
    fvPatchVectorField::operator=(refValue_*patch().nf());
}


SRFSurfaceNormalVelocityFvPatchVectorField::
SRFSurfaceNormalVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    relative_(dict.lookup("relative")),
    refValue_("refValue", dict, p.size())
{
    fvPatchVectorField::operator=(refValue_*patch().nf());
}


SRFSurfaceNormalVelocityFvPatchVectorField::
SRFSurfaceNormalVelocityFvPatchVectorField
(
    const SRFSurfaceNormalVelocityFvPatchVectorField& srfvpvf
)
:
    fixedValueFvPatchVectorField(srfvpvf),
    relative_(srfvpvf.relative_),
    refValue_(srfvpvf.refValue_)
{}


SRFSurfaceNormalVelocityFvPatchVectorField::
SRFSurfaceNormalVelocityFvPatchVectorField
(
    const SRFSurfaceNormalVelocityFvPatchVectorField& srfvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(srfvpvf, iF),
    relative_(srfvpvf.relative_),
    refValue_(srfvpvf.refValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SRFSurfaceNormalVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    vectorField::autoMap(m);
    refValue_.autoMap(m);
}


void SRFSurfaceNormalVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const SRFSurfaceNormalVelocityFvPatchVectorField& tiptf =
        refCast<const SRFSurfaceNormalVelocityFvPatchVectorField>(ptf);

    refValue_.rmap(tiptf.refValue_, addr);
}


void SRFSurfaceNormalVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // If relative, include the effect of the SRF
    if (relative_)
    {
        // Get reference to the SRF model
        const SRF::SRFModel& srf =
            db().lookupObject<SRF::SRFModel>("SRFProperties");

        // Determine patch velocity due to SRF
        const vectorField SRFSurfaceNormalVelocity =
            srf.velocity(patch().Cf());

        operator==(-SRFSurfaceNormalVelocity + refValue_*patch().nf());
    }
    // If absolute, simply supply the inlet value as a fixed value
    else
    {
        operator==(refValue_*patch().nf());
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void SRFSurfaceNormalVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("relative") << relative_ << token::END_STATEMENT << nl;
    refValue_.writeEntry("refValue", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    SRFSurfaceNormalVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
