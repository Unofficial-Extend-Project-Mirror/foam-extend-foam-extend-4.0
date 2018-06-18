/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "chtRcTemperatureFvPatchScalarField.H"
#include "chtRegionCoupleBase.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(p, iF),
    kName_("none"),
    radiation_(false)
{}


Foam::chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
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


Foam::chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
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


Foam::chtRcTemperatureFvPatchScalarField::chtRcTemperatureFvPatchScalarField
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

const Foam::chtRcTemperatureFvPatchScalarField&
Foam::chtRcTemperatureFvPatchScalarField::shadowPatchField() const
{
    if
    (
        !isA<chtRcTemperatureFvPatchScalarField>
        (
            regionCouplingFvPatchScalarField::shadowPatchField()
        )
    )
    {
        FatalErrorIn
        (
            "chtRcTemperatureFvPatchScalarField::shadowPatchField() const"
        )   << "Incorrect shadow patch type for patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " Should be chtRcTemperatureFvPatchScalarField.  Actual type "
            << regionCouplingFvPatchScalarField::shadowPatchField().type()
            << abort(FatalError);
    }

    return dynamic_cast
    <
        const chtRcTemperatureFvPatchScalarField&
    >
    (
        regionCouplingFvPatchScalarField::shadowPatchField()
    );
}


// Return neighbour field given internal cell data
Foam::tmp<Foam::scalarField>
Foam::chtRcTemperatureFvPatchScalarField::patchNeighbourField() const
{
    return regionCouplingFvPatchScalarField::patchNeighbourField
    (
        remoteFieldName()
    );
}


Foam::tmp<Foam::scalarField>
Foam::chtRcTemperatureFvPatchScalarField::Tw() const
{
    return *this;
}


Foam::tmp<Foam::scalarField>
Foam::chtRcTemperatureFvPatchScalarField::Tc() const
{
    return patchInternalField();
}


void Foam::chtRcTemperatureFvPatchScalarField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    const fvPatchScalarField& kpf =
        lookupPatchField<volScalarField, scalar>(kName_);

    if (!isA<chtRegionCoupleBase>(kpf))
    {
        FatalErrorIn
        (
            "void chtRcTemperatureFvPatchScalarField::initEvaluate\n"
            "(\n"
            "    const Pstream::commsTypes commsType\n"
            ")"
        )   << "Incorrect shadow patch type for patch " << this->patch().name()
            << " of field " << dimensionedInternalField().name()
            << ".  Should be derived from chtRegionCoupleBase.  Actual type "
            << kpf.type()
            << abort(FatalError);
    }

    const chtRegionCoupleBase& K = dynamic_cast<const chtRegionCoupleBase&>(kpf);

    *this == K.calcTemperature(*this, shadowPatchField(), K);
}


void Foam::chtRcTemperatureFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchScalarField::evaluate();
}


void Foam::chtRcTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const chtRegionCoupleBase& K =
        refCast<const chtRegionCoupleBase>
        (
            lookupPatchField<volScalarField, scalar>(kName_)
        );

    *this == K.calcTemperature(*this, shadowPatchField(), K);

    fvPatchScalarField::updateCoeffs();
}


Foam::tmp<Foam::scalarField>
Foam::chtRcTemperatureFvPatchScalarField::source() const
{
    const fvPatch& p = patch();

    const scalarField TcOwn = Tc();
    const scalarField TcNei =
        regionCouplePatch().interpolate(shadowPatchField().Tc());
    const scalarField Tw = this->Tw();

    const chtRegionCoupleBase& K =
        dynamic_cast<const chtRegionCoupleBase&>
        (
            p.lookupPatchField<volScalarField, scalar>(kName_)
        );

    const scalarField k = K.kw()*p.deltaCoeffs();
    const scalarField kOwn = K.korig();

    return kOwn*(Tw - TcOwn) - k*(TcNei - TcOwn);
}


void Foam::chtRcTemperatureFvPatchScalarField::manipulateMatrix
(
    fvScalarMatrix& matrix
)
{
    const fvPatch& p = patch();
    const scalarField& magSf = p.magSf();
    const labelList& cellLabels = p.faceCells();
    scalarField& source = matrix.source();

    scalarField s = this->source();

    forAll(cellLabels, i)
    {
        source[cellLabels[i]] += s[i]*magSf[i];
    }
}


void Foam::chtRcTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("K") << kName_ << token::END_STATEMENT << nl;
    os.writeKeyword("radiation") << radiation_ << token::END_STATEMENT << nl;
    os.writeKeyword("remoteField")
        << remoteFieldName() << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    chtRcTemperatureFvPatchScalarField
);

} // End namespace Foam


// ************************************************************************* //
