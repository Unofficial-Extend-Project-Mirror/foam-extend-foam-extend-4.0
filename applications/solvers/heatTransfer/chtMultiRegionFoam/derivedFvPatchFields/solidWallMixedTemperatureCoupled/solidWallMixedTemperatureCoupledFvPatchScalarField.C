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

#include "solidWallMixedTemperatureCoupledFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::
solidWallMixedTemperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    coupleManager_(p),
    KName_("undefined-K")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
    this->fixesValue_ = true;
}


Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::
solidWallMixedTemperatureCoupledFvPatchScalarField
(
    const solidWallMixedTemperatureCoupledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    coupleManager_(ptf.coupleManager_),
    KName_(ptf.KName_),
    fixesValue_(ptf.fixesValue_)
{}


Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::
solidWallMixedTemperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    coupleManager_(p, dict),
    KName_(dict.lookup("K"))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    refValue() = static_cast<scalarField>(*this);
    refGrad() = 0.0;
    valueFraction() = 1.0;
    fixesValue_ = true;
}


Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::
solidWallMixedTemperatureCoupledFvPatchScalarField
(
    const solidWallMixedTemperatureCoupledFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    coupleManager_(wtcsf.coupleManager_),
    KName_(wtcsf.KName_),
    fixesValue_(wtcsf.fixesValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvPatchScalarField&
Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::K() const
{
    return this->patch().lookupPatchField<volScalarField, scalar>(KName_);
}


void Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    tmp<scalarField> intFld = patchInternalField();

    label nFixed = 0;

    // Like snGrad but bypass switching on refValue/refGrad.
    tmp<scalarField> normalGradient =
        (*this-intFld())
      * patch().deltaCoeffs();

    if (debug)
    {
        Info<< "solidWallMixedTemperatureCoupledFvPatchScalarField::"
            << "updateCoeffs() :"
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    forAll(*this, i)
    {
        // if outgoing flux use fixed value.
        if (intFld()[i] > operator[](i))
        {
            this->refValue()[i] = operator[](i);
            this->refGrad()[i] = 0.0;   // not used
            this->valueFraction()[i] = 1.0;
            nFixed++;
        }
        else
        {
            this->refValue()[i] = 0.0;  // not used
            this->refGrad()[i] = normalGradient()[i];
            this->valueFraction()[i] = 0.0;
        }
    }

    reduce(nFixed, sumOp<label>());

    fixesValue_ = (nFixed > 0);

    if (debug)
    {
        label nTotSize = returnReduce(this->size(), sumOp<label>());

        Info<< "solidWallMixedTemperatureCoupledFvPatchScalarField::"
            << "updateCoeffs() : Out of " << nTotSize
            << " fixedBC:" << nFixed
            << " gradient:" << nTotSize-nFixed << endl;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if (!coupleManager_.regionOwner())
    {
        // I am the last one to evaluate.

        tmp<scalarField> intFld = patchInternalField();

        const fvPatch& nbrPatch = coupleManager_.neighbourPatch();

        solidWallMixedTemperatureCoupledFvPatchScalarField& nbrField =
        refCast<solidWallMixedTemperatureCoupledFvPatchScalarField>
        (
            const_cast<fvPatchField<scalar>&>
            (
                coupleManager_.neighbourPatchField<scalar>()
            )
        );
        tmp<scalarField> nbrIntFld = nbrField.patchInternalField();
        tmp<scalarField> myKDelta = K()*patch().deltaCoeffs();
        tmp<scalarField> nbrKDelta = nbrField.K()*nbrPatch.deltaCoeffs();

        // Calculate common wall temperature and assign to both sides
        scalarField::operator=
        (
            (myKDelta()*intFld + nbrKDelta()*nbrIntFld)
          / (myKDelta() + nbrKDelta())
        );

        nbrField.scalarField::operator=(*this);

        if (debug)
        {
            Info<< "Setting master and slave to wall temperature "
                << " min:" << gMin(*this)
                << " max:" << gMax(*this)
                << " avg:" << gAverage(*this)
                << endl;
        }
    }

    fvPatchScalarField::evaluate();
}


void Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    coupleManager_.writeEntries(os);
    os.writeKeyword("K") << KName_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    solidWallMixedTemperatureCoupledFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
