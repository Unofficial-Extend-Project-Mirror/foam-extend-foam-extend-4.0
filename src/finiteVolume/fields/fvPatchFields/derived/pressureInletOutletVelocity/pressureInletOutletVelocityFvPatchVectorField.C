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

#include "pressureInletOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pressureInletOutletVelocityFvPatchVectorField::
pressureInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    phiName_("phi")
{
    refValue() = vector::zero;
    refGrad() = vector::zero;
    valueFraction() = symmTensor::zero;
}


pressureInletOutletVelocityFvPatchVectorField::
pressureInletOutletVelocityFvPatchVectorField
(
    const pressureInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
{
    if (ptf.tangentialVelocity_.size())
    {
        tangentialVelocity_ = mapper(ptf.tangentialVelocity_);
    }
}


pressureInletOutletVelocityFvPatchVectorField::
pressureInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    if (dict.found("tangentialVelocity"))
    {
        setTangentialVelocity
        (
            vectorField("tangentialVelocity", dict, p.size())
        );
    }
    else
    {
        refValue() = vector::zero;
    }

    refGrad() = vector::zero;
    valueFraction() = symmTensor::zero;
}


pressureInletOutletVelocityFvPatchVectorField::
pressureInletOutletVelocityFvPatchVectorField
(
    const pressureInletOutletVelocityFvPatchVectorField& pivpvf
)
:
    directionMixedFvPatchVectorField(pivpvf),
    phiName_(pivpvf.phiName_),
    tangentialVelocity_(pivpvf.tangentialVelocity_)
{}


pressureInletOutletVelocityFvPatchVectorField::
pressureInletOutletVelocityFvPatchVectorField
(
    const pressureInletOutletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    tangentialVelocity_(pivpvf.tangentialVelocity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pressureInletOutletVelocityFvPatchVectorField::
setTangentialVelocity(const vectorField& tangentialVelocity)
{
    tangentialVelocity_ = tangentialVelocity;
    vectorField n = patch().nf();
    refValue() = tangentialVelocity_ - n*(n & tangentialVelocity_);
}


void pressureInletOutletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
    if (tangentialVelocity_.size())
    {
        tangentialVelocity_.autoMap(m);
    }
}


void pressureInletOutletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);

    if (tangentialVelocity_.size())
    {
        const pressureInletOutletVelocityFvPatchVectorField& tiptf =
            refCast<const pressureInletOutletVelocityFvPatchVectorField>(ptf);

        tangentialVelocity_.rmap(tiptf.tangentialVelocity_, addr);
    }
}


void pressureInletOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (this->db().objectRegistry::found(phiName_))
    {
        const fvsPatchField<scalar>& phip =
            lookupPatchField<surfaceScalarField, scalar>(phiName_);

        valueFraction() = neg(phip)*(I - sqr(patch().nf()));
    }
    else
    {
        InfoIn
        (
            "pressureInletOutletVelocityFvPatchVectorField::updateCoeffs()"
        )   << "Cannot find phi.  Return" << endl;

        valueFraction() = symmTensor::one;
    }

    directionMixedFvPatchVectorField::updateCoeffs();
}


void pressureInletOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    if (tangentialVelocity_.size())
    {
        tangentialVelocity_.writeEntry("tangentialVelocity", os);
    }
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void pressureInletOutletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    vectorField normalValue = transform(valueFraction(), refValue());
    vectorField transformGradValue = transform(I - valueFraction(), pvf);
    fvPatchField<vector>::operator=(normalValue + transformGradValue);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    pressureInletOutletVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
