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

Description

\*---------------------------------------------------------------------------*/

#include "cohesiveLawFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "rheologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cohesiveLawFvPatchVectorField::cohesiveLawFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    cohesiveLawPtr_(NULL),
    relaxationFactor_(1.0),
    traction_(p.size(), vector::zero)
{}


cohesiveLawFvPatchVectorField::cohesiveLawFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    cohesiveLawPtr_
    (
        cohesiveLaw::New(dict.lookup("cohesiveLaw"), dict).ptr()
    ),
    relaxationFactor_(readScalar(dict.lookup("relaxationFactor"))),
    traction_(p.size(), vector::zero)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


cohesiveLawFvPatchVectorField::cohesiveLawFvPatchVectorField
(
    const cohesiveLawFvPatchVectorField& cpf
)
:
    fixedGradientFvPatchVectorField(cpf),
    cohesiveLawPtr_(cpf.cohesiveLawPtr_->clone().ptr()),
    relaxationFactor_(cpf.relaxationFactor_),
    traction_(cpf.traction_)
{}


cohesiveLawFvPatchVectorField::cohesiveLawFvPatchVectorField
(
    const cohesiveLawFvPatchVectorField& cpf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(cpf, p, iF, mapper),
    cohesiveLawPtr_(cpf.cohesiveLawPtr_->clone().ptr()),
    relaxationFactor_(cpf.relaxationFactor_),
    traction_(cpf.traction_, mapper)
{}


cohesiveLawFvPatchVectorField::cohesiveLawFvPatchVectorField
(
    const cohesiveLawFvPatchVectorField& cpf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(cpf, iF),
    cohesiveLawPtr_(cpf.cohesiveLawPtr_->clone().ptr()),
    relaxationFactor_(cpf.relaxationFactor_),
    traction_(cpf.traction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const cohesiveLaw& cohesiveLawFvPatchVectorField::law() const
{
    if (!cohesiveLawPtr_)
    {
        FatalErrorIn
        (
            "const cohesiveLaw& cohesiveLawFvPatchVectorField::law() const"
        )   << "Law pointer not set" << abort(FatalError);
    }

    return *cohesiveLawPtr_;
}


void cohesiveLawFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (cohesiveLawPtr_ == NULL)
    {
        FatalErrorIn("cohesiveFvPatchVectorField::autoMap")
            << "NULL cohesive law"
            << abort(FatalError);
    }

    fixedGradientFvPatchVectorField::autoMap(m);

    traction_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void cohesiveLawFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const cohesiveLawFvPatchVectorField& dmptf =
        refCast<const cohesiveLawFvPatchVectorField>(ptf);

    // No need to grab the cohesive zone pointer more than once
    if (!cohesiveLawPtr_)
    {
        cohesiveLawPtr_ = dmptf.cohesiveLawPtr_->clone().ptr();

        relaxationFactor_ = dmptf.relaxationFactor_;
    }

    traction_.rmap(dmptf.traction_, addr);
}


// Update the coefficients associated with the patch field
void cohesiveLawFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Looking up rheology

    const fvPatchField<scalar>& mu =
        patch().lookupPatchField<volScalarField, scalar>("mu");

    const fvPatchField<scalar>& lambda =
      patch().lookupPatchField<volScalarField, scalar>("lambda");

    vectorField n = patch().nf();

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>("grad(U)");

    // Patch displacement
    const vectorField& U = *this;

    // Patch stress
    tensorField sigma = mu*(gradU + gradU.T()) + I*(lambda*tr(gradU));

    // Normal stress component
    scalarField sigmaN = (n & (n & sigma));

    scalarField delta = -(n & U);

    label sizeByTwo = patch().size()/2;

    for(label i = 0; i < sizeByTwo; i++)
    {
        scalar tmp = delta[i];
        delta[i] += delta[sizeByTwo + i];
        delta[sizeByTwo + i] += tmp;
    }

    forAll (traction_, faceI)
    {
        if (delta[faceI] < 0)
        {
            // Return from traction to symmetryPlane??
            traction_[faceI] = law().sigmaMax().value()*n[faceI];
        }
        else if(delta[faceI] > law().deltaC().value())
        {
            // Traction free
            traction_[faceI] = vector::zero;
        }
        else
        {
            // Calculate cohesive traction from cohesive zone model
            traction_[faceI] = law().traction(delta[faceI])*n[faceI];
        }
    }

    gradient() =
    (
        traction_
      - (n & (mu*gradU.T() - (mu + lambda)*gradU))
      - n*lambda*tr(gradU)
    )/(2.0*mu + lambda);

    fixedGradientFvPatchVectorField::updateCoeffs();
}


// Write
void cohesiveLawFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    os.writeKeyword("cohesiveLaw") << law().type() 
        << token::END_STATEMENT << nl;
    os.writeKeyword("relaxationFactor") << relaxationFactor_
        << token::END_STATEMENT << nl;
    law().writeDict(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, cohesiveLawFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
