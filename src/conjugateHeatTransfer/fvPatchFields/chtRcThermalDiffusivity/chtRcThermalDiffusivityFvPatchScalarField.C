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
#include "basicThermo.H"

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

    if
    (
        dimensionedInternalField().dimensions()
     == dimensionSet(1, -1, -1, 0, 0, 0, 0)
    )
    {
        const label patchi = patch().index();

        const basicThermo& thermo = db().lookupObject<basicThermo>
        (
            "thermophysicalProperties"
        );

        const scalarField Tw = lookupPatchField<volScalarField, scalar>("T");

        // Note: BUG work-around
        // Selection of enthaly field cannot just rely on dimensions
        // Appropriate function should be added into basicThermo
        // which identifies the active form or enthalpy at run-time
        // in order to avoid nonImplemented virtual function call
        // Fix it.  HJ, 15/Jun/2018
        if (this->db().objectRegistry::foundObject<volScalarField>("h"))
        {
            const chtRcTemperatureFvPatchScalarField& h =
                refCast<const chtRcTemperatureFvPatchScalarField>
                (
                    thermo.h().boundaryField()[patchi]
                );

            *this == calcThermalDiffusivity(*this, shadowPatchField(), h)
            /thermo.Cp(Tw, patchi);
        }
        else if
        (
            this->db().objectRegistry::foundObject<volScalarField>("hs")
        )
        {
            const chtRcTemperatureFvPatchScalarField& hs =
                refCast<const chtRcTemperatureFvPatchScalarField>
                (
                    thermo.hs().boundaryField()[patchi]
                );

            *this == calcThermalDiffusivity(*this, shadowPatchField(), hs)
                /thermo.Cp(Tw, patchi);
        }
        else
        {
            FatalErrorIn
            (
                "void fixedFluxPressureFvPatchScalarField::updateCoeffs()"
            )   << "Cannot find enthalpy for field "
                << dimensionedInternalField().name() << " on patch "
                << patch().name()
                << abort(FatalError);
        }
    }
}


Foam::tmp<Foam::scalarField>
Foam::chtRcThermalDiffusivityFvPatchScalarField::calcThermalDiffusivity
(
    const chtRegionCoupleBase& owner,
    const chtRegionCoupleBase& neighbour,
    const chtRcTemperatureFvPatchScalarField& TwOwn
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

    const scalarField fOwn = owner.forig();
    const scalarField TcOwn = TwOwn.Tc();

    scalarField fNei(p.size());
    scalarField TcNei(p.size());

    scalarField Qr(p.size(), 0.0);
    scalarField fourQro(p.size(), 0.0);

    if (TwOwn.radiation())
    {
        Qr += p.lookupPatchField<volScalarField, scalar>("Qr");
        fourQro += 4.0*radiation::sigmaSB.value()*pow4(TwOwn.Tw());
    }

    {
        Field<VectorN<scalar, 4> > lData(neighbour.size());

        const scalarField lfNei = neighbour.forig();
        const scalarField lTcNei = TwOwn.shadowPatchField().Tc();

        forAll(lData, facei)
        {
            lData[facei][0] = lTcNei[facei];
            lData[facei][1] = lfNei[facei];
        }

        if(TwOwn.shadowPatchField().radiation())
        {
            const scalarField& lTwNei = TwOwn.shadowPatchField().Tw();
            const scalarField& lQrNei =
                owner.lookupShadowPatchField<volScalarField, scalar>("Qr");

            forAll (lData, facei)
            {
                lData[facei][2] = lTwNei[facei];
                lData[facei][3] = lQrNei[facei];
            }
        }
        else
        {
            forAll (lData, facei)
            {
                lData[facei][2] = 0.0;
                lData[facei][3] = 0.0;
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

    tmp<scalarField> kTmp(new scalarField(p.size()));
    scalarField& k = kTmp();

    k = kOwn*(kNei + Qr/stabilise(TcNei - TcOwn, SMALL));
    k /= p.deltaCoeffs()*(kOwn + kNei);

    forAll (k, facei)
    {
        k[facei] = max(min(k[facei], 100*kHarm[facei]), 0.01*kHarm[facei]);
    }

    return kTmp;
}


Foam::tmp<Foam::scalarField>
Foam::chtRcThermalDiffusivityFvPatchScalarField::calcTemperature
(
    const chtRcTemperatureFvPatchScalarField& TwOwn,
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

    const scalarField fOwn = ownerK.forig();
    const scalarField TcOwn = TwOwn.Tc();

    scalarField fNei(p.size());
    scalarField TcNei(p.size());

    scalarField Qr(p.size(), 0.0);
    scalarField fourQro(p.size(), 0.0);

    if (TwOwn.radiation())
    {
        Qr += p.lookupPatchField<volScalarField, scalar>("Qr");
        fourQro += 4.0*radiation::sigmaSB.value()*pow4(TwOwn.Tw());
    }

    {
        Field<VectorN<scalar, 4> > lData(neighbour.size());

        const scalarField lfNei = ownerK.shadowPatchField().forig();
        const scalarField lTcNei = TwOwn.shadowPatchField().Tc();

        forAll (lData, facei)
        {
            lData[facei][0] = lTcNei[facei];
            lData[facei][1] = lfNei[facei];
        }

        if (TwOwn.shadowPatchField().radiation())
        {
            const scalarField& lTwNei = TwOwn.shadowPatchField().Tw();
            const scalarField& lQrNei =
                TwOwn.lookupShadowPatchField<volScalarField, scalar>("Qr");

            forAll (lData, facei)
            {
                lData[facei][2] = lTwNei[facei];
                lData[facei][3] = lQrNei[facei];
            }
        }
        else
        {
            forAll (lData, facei)
            {
                lData[facei][2] = 0.0;
                lData[facei][3] = 0.0;
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

    // Do interpolation
    harmonic<scalar> interp(mesh);
    scalarField weights = interp.weights(fOwn, fNei, p);
    const scalarField kHarm = weights*fOwn + (1.0 - weights)*fNei;

    const scalarField kOwn = fOwn/(1.0 - p.weights())/mld.magDelta(p.index());
    const scalarField kNei = fNei/p.weights()/mld.magDelta(p.index());

    tmp<scalarField> TwTmp(new scalarField(TwOwn.Tw()));
    scalarField& Tw = TwTmp();

    Tw = (Qr + kOwn*TcOwn + kNei*TcNei)/(kOwn + kNei);

    scalarField q1 = (Tw - TcOwn)*kOwn;

    scalarField q2 = (TcNei - Tw)*kNei;

    scalarField q3 = (TcNei - TcOwn)*ownerK*p.deltaCoeffs();

    return TwTmp;
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
