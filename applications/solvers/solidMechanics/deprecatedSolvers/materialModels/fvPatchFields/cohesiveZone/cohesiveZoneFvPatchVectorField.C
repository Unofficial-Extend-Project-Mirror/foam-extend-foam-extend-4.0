/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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

Class
    cohesiveZoneFvPatchVectorField

Description

\*---------------------------------------------------------------------------*/

#include "cohesiveZoneFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "rheologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cohesiveZoneFvPatchVectorField::cohesiveZoneFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    UName_("undefined"),
    rheologyName_("undefined"),
    cohesiveLawPtr_(NULL),
    relaxationFactor_(1.0)
{}


cohesiveZoneFvPatchVectorField::cohesiveZoneFvPatchVectorField
(
    const cohesiveZoneFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    rheologyName_(ptf.rheologyName_),
    cohesiveLawPtr_(ptf.cohesiveLawPtr_),
    relaxationFactor_(ptf.relaxationFactor_)
{}


cohesiveZoneFvPatchVectorField::cohesiveZoneFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    UName_(dict.lookup("U")),
    rheologyName_(dict.lookup("rheology")),
    cohesiveLawPtr_
    (
        cohesiveLaw::New(dict.lookup("cohesiveLaw"), dict).ptr()
    ),
    relaxationFactor_(readScalar(dict.lookup("relaxationFactor")))
{
    if (dict.found("refValue"))
    {
        this->refValue() = vectorField("refValue", dict, p.size());
    }
    else
    {
        this->refValue() = vector::zero;
    }
    
    if (dict.found("refGradient"))
    {
        this->refGrad() = vectorField("refGradient", dict, p.size());
    }
    else
    {
        this->refGrad() = vector::zero;
    }

    if (dict.found("valueFraction"))
    {
        this->valueFraction() = 
            symmTensorField("valueFraction", dict, p.size());
    }
    else
    {
        vectorField n = patch().nf();

        this->valueFraction() = sqr(n);
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        Field<vector> normalValue = transform(valueFraction(), refValue());

        Field<vector> gradValue =
            this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

        Field<vector> transformGradValue =
            transform(I - valueFraction(), gradValue);

        Field<vector>::operator=(normalValue + transformGradValue);
    }
}


cohesiveZoneFvPatchVectorField::cohesiveZoneFvPatchVectorField
(
    const cohesiveZoneFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    UName_(ptf.UName_),
    rheologyName_(ptf.rheologyName_),
    cohesiveLawPtr_(ptf.cohesiveLawPtr_),
    relaxationFactor_(ptf.relaxationFactor_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void cohesiveZoneFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (cohesiveLawPtr_ == NULL)
    {
        FatalErrorIn("cohesiveZoneFvPatchVectorField::autoMap")
            << "NULL cohesive law"
            << abort(FatalError);
    }

    directionMixedFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void cohesiveZoneFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);

    const cohesiveZoneFvPatchVectorField& dmptf =
        refCast<const cohesiveZoneFvPatchVectorField>(ptf);

    // No need to grab the cohesive zone pointer more than once
    if (!cohesiveLawPtr_)
    {
        cohesiveLawPtr_ = dmptf.cohesiveLawPtr_->clone().ptr();

        relaxationFactor_ = dmptf.relaxationFactor_;
    }
}


void cohesiveZoneFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Looking up rheology
    const rheologyModel& rheology =
        this->db().objectRegistry::lookupObject<rheologyModel>(rheologyName_);

    const scalarField mu = 
        rheology.mu()().boundaryField()[patch().index()];

    const scalarField lambda =
        rheology.lambda()().boundaryField()[patch().index()];

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" +UName_ + ")"
        );

    // Patch displacement
    const vectorField& U = *this;

    // Patch stress
    tensorField sigma = mu*(gradU + gradU.T()) + I*(lambda*tr(gradU));

    // Patch normal
    vectorField n = patch().nf();

    // Normal stress component
    scalarField sigmaN = (n&(n&sigma));

    // Chech crack propagation
    forAll(sigmaN, faceI)
    {
        vector cohesiveTraction = vector::zero;

        if
        (
            (magSqr(valueFraction()[faceI]) > 1-SMALL)
         && (sigmaN[faceI] >= law().sigmaMax().value())
        )
        {
            // Switch to full traction boundary condition
            valueFraction()[faceI] = symmTensor::zero;

            Info << "Crack started at face: " << faceI << endl;

            // Cohesive traction
            cohesiveTraction = n[faceI]*law().sigmaMax().value();
        }
        else if(magSqr(valueFraction()[faceI]) < SMALL)
        {
            // Normal displacement
            scalar Un = -(n[faceI]&U[faceI]);

            if(Un < 0)
            {
                // Return from traction to symmetryPlane
                refValue()[faceI] = vector::zero;
                refGrad() = vector::zero;
                valueFraction()[faceI] = sqr(n[faceI]);
                Info << "Face removed from crack: " << faceI << endl;
            }
            else if(Un > law().deltaC().value()/2)
            {
                // Traction free
                cohesiveTraction = vector::zero;
            }
            else
            {
                // Calculate cohesive traction from cohesive zone model
                cohesiveTraction = law().traction(2*Un)*n[faceI];
            }
        }

        if(magSqr(valueFraction()[faceI]) < SMALL)
        {
            cohesiveTraction = 
                relaxationFactor_*cohesiveTraction 
              + (1.0 - relaxationFactor_)*sigmaN[faceI]*n[faceI];

            refGrad()[faceI] =
            (
                cohesiveTraction
              - (
                    n[faceI] 
                  & (
                        mu[faceI]*gradU[faceI].T() 
                      - (mu[faceI] + lambda[faceI])*gradU[faceI]
                    )
                )
              - n[faceI]*lambda[faceI]*tr(gradU[faceI])
            )
           /(2.0*mu[faceI] + lambda[faceI]);
        }
    }

    directionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void cohesiveZoneFvPatchVectorField::write(Ostream& os) const
{
    directionMixedFvPatchVectorField::write(os);
    os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rheology") << rheologyName_ << token::END_STATEMENT << nl;
    os.writeKeyword("cohesiveLaw") << law().type() 
        << token::END_STATEMENT << nl;
    os.writeKeyword("relaxationFactor") << relaxationFactor_
        << token::END_STATEMENT << nl;
    law().writeDict(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, cohesiveZoneFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
