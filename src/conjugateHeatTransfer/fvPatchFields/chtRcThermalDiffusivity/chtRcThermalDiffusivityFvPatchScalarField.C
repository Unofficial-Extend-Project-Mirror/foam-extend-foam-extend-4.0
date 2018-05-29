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

#include "chtRcThermalDiffusivityFvPatchScalarField.H"
#include "chtRcThermalDiffusivitySlaveFvPatchScalarField.H"
#include "chtRcTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "harmonic.H"
#include "radiationConstants.H"
#include "VectorN.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chtRcThermalDiffusivityFvPatchScalarField::
chtRcThermalDiffusivityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRegionCoupleBase(p, iF)
{}


Foam::chtRcThermalDiffusivityFvPatchScalarField::
chtRcThermalDiffusivityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    chtRegionCoupleBase(p, iF, dict)
{}


Foam::chtRcThermalDiffusivityFvPatchScalarField::
chtRcThermalDiffusivityFvPatchScalarField
(
    const chtRcThermalDiffusivityFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    chtRegionCoupleBase(ptf, p, iF, mapper)
{}


Foam::chtRcThermalDiffusivityFvPatchScalarField::
chtRcThermalDiffusivityFvPatchScalarField
(
    const chtRcThermalDiffusivityFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRegionCoupleBase(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::chtRcThermalDiffusivityFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchScalarField::evaluate();
}


void Foam::chtRcThermalDiffusivityFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    calcThermalDiffusivity(*this, shadowPatchField());
}


void
Foam::chtRcThermalDiffusivityFvPatchScalarField::calcThermalDiffusivity
(
    chtRegionCoupleBase& owner,
    const chtRegionCoupleBase& neighbour
) const
{
    if (debug)
    {
        InfoIn
        (
            "chtRcThermalDiffusivityFvPatchScalarField::calcThermalDiffusivity"
        )   << "for field " << this->dimensionedInternalField().name()
            << " in " << this->patch().boundaryMesh().mesh().name()
            << endl;
    }

    const fvPatch& p = owner.patch();
    const fvMesh& mesh = p.boundaryMesh().mesh();
    const magLongDelta& mld = magLongDelta::New(mesh);

    const chtRcTemperatureFvPatchScalarField& TwOwn =
        dynamic_cast<const chtRcTemperatureFvPatchScalarField&>
        (
            p.lookupPatchField<volScalarField, scalar>("T")
        );
    scalarField& k = owner;
    const scalarField& fOwn = owner.originalPatchField();
    const scalarField TcOwn = TwOwn.patchInternalField();

    scalarField fNei(p.size());
    scalarField TcNei(p.size());

    scalarField Qr(p.size(), 0.0);
    scalarField fourQro(p.size(), 0.0);

    if (TwOwn.radiation())
    {
        Qr += p.lookupPatchField<volScalarField, scalar>("Qr");
        fourQro += 4.0*radiation::sigmaSB.value()*pow4(TwOwn);
    }

    {
        Field<VectorN<scalar, 4> > lData
        (
            neighbour.size(),
            pTraits<VectorN<scalar, 4> >::zero
        );

        const scalarField& lfNei = neighbour.originalPatchField();
        scalarField lTcNei = TwOwn.shadowPatchField().patchInternalField();

        forAll (lData, facei)
        {
            lData[facei][0] = lTcNei[facei];
            lData[facei][1] = lfNei[facei];
        }

        if (TwOwn.shadowPatchField().radiation())
        {
            const scalarField& lQrNei =
                owner.lookupShadowPatchField<volScalarField, scalar>("Qr");
            const scalarField& lTwNei = TwOwn.shadowPatchField();

            forAll (lData, facei)
            {
                lData[facei][2] = lTwNei[facei];
                lData[facei][3] = lQrNei[facei];
            }
        }

        const Field<VectorN<scalar, 4> > iData =
            owner.regionCouplePatch().interpolate(lData);

        forAll (iData, facei)
        {
            TcNei[facei] = iData[facei][0];
            fNei[facei] = iData[facei][1];
        }

        if (TwOwn.shadowPatchField().radiation())
        {
            forAll (iData, facei)
            {
                Qr[facei] += iData[facei][3];
                fourQro[facei] +=
                    4.0*radiation::sigmaSB.value()*pow4(iData[facei][2]);
            }
        }
    }

    // Do interpolation
    harmonic<scalar> interp(mesh);
    const scalarField weights = interp.weights(fOwn, fNei, p);
    const scalarField kHarm = weights*fOwn + (1.0 - weights)*fNei;

    const scalarField kOwn = fOwn/(1.0 - p.weights())/mld.magDelta(p.index());
    const scalarField kNei = fNei/p.weights()/mld.magDelta(p.index());

    k = kOwn*(TwOwn*(kNei*(TcNei - TcOwn) + Qr + fourQro) - TcOwn*fourQro);
    k /= stabilise((fourQro + TwOwn*(kOwn + kNei))*(TcNei - TcOwn), SMALL);
    k /= p.deltaCoeffs();

    //Info << "k = " << k << endl;

    forAll (k, facei)
    {
        k[facei] = max(min(k[facei], 100*kHarm[facei]), 0.01*kHarm[facei]);
    }

    owner.fvPatchScalarField::updateCoeffs();
}


void
Foam::chtRcThermalDiffusivityFvPatchScalarField::calcTemperature
(
    chtRcTemperatureFvPatchScalarField& TwOwn,
    const chtRcTemperatureFvPatchScalarField& neighbour,
    const chtRegionCoupleBase& ownerK
) const
{
    if (debug)
    {
        InfoIn
        (
            "chtRcThermalDiffusivityFvPatchScalarField::calcTemperature"
        )   << "for field " << this->dimensionedInternalField().name()
            << " in " << this->patch().boundaryMesh().mesh().name()
            << endl;
    }

    const fvPatch& p = TwOwn.patch();
    const fvMesh& mesh = p.boundaryMesh().mesh();
    const magLongDelta& mld = magLongDelta::New(mesh);

    const scalarField& fOwn = ownerK.originalPatchField();
    const scalarField TcOwn = TwOwn.patchInternalField();

    scalarField fNei(p.size());
    scalarField TcNei(p.size());

    scalarField Qr(p.size(), 0.0);
    scalarField fourQro(p.size(), 0.0);

    if (TwOwn.radiation())
    {
        Qr += p.lookupPatchField<volScalarField, scalar>("Qr");
        fourQro += 4.0*radiation::sigmaSB.value()*pow4(TwOwn);
    }

    {
        Field<VectorN<scalar, 4> > lData
        (
            neighbour.size(),
            pTraits<VectorN<scalar, 4> >::zero
        );

        const scalarField& lfNei =
            ownerK.shadowPatchField().originalPatchField();
        scalarField lTcNei =
            TwOwn.shadowPatchField().patchInternalField();

        forAll (lData, facei)
        {
            lData[facei][0] = lTcNei[facei];
            lData[facei][1] = lfNei[facei];
        }

        if (TwOwn.shadowPatchField().radiation())
        {
            const scalarField& lTwNei = TwOwn.shadowPatchField();
            const scalarField& lQrNei =
                TwOwn.lookupShadowPatchField<volScalarField, scalar>("Qr");

            forAll (lData, facei)
            {
                lData[facei][2] = lTwNei[facei];
                lData[facei][3] = lQrNei[facei];
            }
        }

        const Field<VectorN<scalar, 4> > iData =
            TwOwn.regionCouplePatch().interpolate(lData);

        forAll (iData, facei)
        {
            TcNei[facei] = iData[facei][0];
            fNei[facei] = iData[facei][1];
        }

        if (TwOwn.shadowPatchField().radiation())
        {
            forAll (iData, facei)
            {
                fourQro[facei] +=
                    4.0*radiation::sigmaSB.value()*pow4(iData[facei][2]);
                Qr[facei] += iData[facei][3];
            }
        }
    }

    const scalarField kOwn = fOwn/(1.0 - p.weights())/mld.magDelta(p.index());
    const scalarField kNei = fNei/p.weights()/mld.magDelta(p.index());

    TwOwn *=
        (fourQro + Qr + kOwn*TcOwn + kNei*TcNei)
       /(TwOwn*(kOwn + kNei) + fourQro);

    TwOwn.fvPatchScalarField::updateCoeffs();
}


void Foam::chtRcThermalDiffusivityFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("remoteField")
        << remoteFieldName() << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    makePatchTypeField
    (
        fvPatchScalarField,
        chtRcThermalDiffusivityFvPatchScalarField
    );

} // End namespace Foam


// ************************************************************************* //
