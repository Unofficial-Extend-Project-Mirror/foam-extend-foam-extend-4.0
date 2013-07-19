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

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "chtRcTemperatureFvPatchScalarField.H"
#include "chtRegionCoupleBase.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "magLongDelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(p, iF),
    kName_("none"),
    radiation_(false)
{}


chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    regionCouplingFvPatchScalarField(p, iF, dict),
    kName_(dict.lookup("K")),
    radiation_(readBool(dict.lookup("radiation")))
{}


chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
(
    const chtRcTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    regionCouplingFvPatchScalarField(ptf, p, iF, mapper),
    kName_(ptf.kName_),
    radiation_(ptf.radiation_)
{}


chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
(
    const chtRcTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(ptf, iF),
    kName_(ptf.kName_),
    radiation_(ptf.radiation_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// Return a named shadow patch field
const chtRcTemperatureFvPatchScalarField&
chtRcTemperatureFvPatchScalarField::shadowPatchField() const
{
    return dynamic_cast
    <
        const chtRcTemperatureFvPatchScalarField&
    >
    (
        regionCouplingFvPatchScalarField::shadowPatchField()
    );
}


void chtRcTemperatureFvPatchScalarField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    const chtRegionCoupleBase& K =
        dynamic_cast<const chtRegionCoupleBase&>
        (
            patch().lookupPatchField<volScalarField, scalar>(kName_)
        );

    K.calcTemperature(*this, shadowPatchField(), K);
}


void chtRcTemperatureFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchScalarField::evaluate();
}


void chtRcTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fvPatchScalarField::updateCoeffs();
}


tmp<scalarField> chtRcTemperatureFvPatchScalarField::source() const
{
    const fvPatch& p = patch();
    const magLongDelta& mld = magLongDelta::New(p.boundaryMesh().mesh());

    const scalarField TcOwn = patchInternalField();
    const scalarField TcNei = patchNeighbourField();
    const scalarField& Tw = *this;

    const chtRegionCoupleBase& K =
        dynamic_cast<const chtRegionCoupleBase&>
        (
            p.lookupPatchField<volScalarField, scalar>(kName_)
        );

    const scalarField k = K*p.deltaCoeffs();

    const scalarField kOwn =
        K.originalPatchField()/(1 - p.weights())/mld.magDelta(p.index());

    return kOwn*(Tw - TcOwn) - k*(TcNei - TcOwn);
}


void chtRcTemperatureFvPatchScalarField::manipulateMatrix
(
    fvScalarMatrix& matrix
)
{
    const fvPatch& p = patch();
    const scalarField& magSf = p.magSf();
    const labelList& cellLabels = p.faceCells();
    scalarField& source = matrix.source();

    scalarField s = this->source();

    //Info << "s = " << s << " Sum s = " << sum(s*p.magSf()) << endl;
    //Info << "Sum s = " << sum(s*p.magSf()) << endl;

    forAll(cellLabels, i)
    {
        source[cellLabels[i]] += s[i]*magSf[i];
    }
}


void chtRcTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("K") << kName_ << token::END_STATEMENT << nl;
    os.writeKeyword("radiation") << radiation_ << token::END_STATEMENT << nl;
    os.writeKeyword("remoteField")
        << remoteFieldName() << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    chtRcTemperatureFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
