/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "directMappedPatchBase.H"
#include "regionProperties.H"
#include "basicThermo.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_("undefined-neighbourFieldName"),
    KappaName_("undefined-Kappa")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    neighbourFieldName_(ptf.neighbourFieldName_),
    KappaName_(ptf.KappaName_)
{}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_(dict.lookup("neighbourFieldName")),
    KappaName_(dict.lookup("Kappa"))
{
    if (!isA<directMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::"
            "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << directMappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    neighbourFieldName_(wtcsf.neighbourFieldName_),
    KappaName_(wtcsf.KappaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField>
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::Kappa() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (KappaName_ == "none")
    {
        const compressible::RASModel& model =
            db().lookupObject<compressible::RASModel>("RASProperties");

        tmp<volScalarField> talpha = model.alphaEff();

        const basicThermo& thermo =
            db().lookupObject<basicThermo>("thermophysicalProperties");

        return
            talpha().boundaryField()[patch().index()]
           *thermo.Cp()().boundaryField()[patch().index()];
    }
    else if (mesh.objectRegistry::foundObject<volScalarField>(KappaName_))
    {
        return lookupPatchField<volScalarField, scalar>(KappaName_);
    }
    else if (mesh.objectRegistry::foundObject<volSymmTensorField>(KappaName_))
    {
        const symmTensorField& KappaWall =
            lookupPatchField<volSymmTensorField, scalar>(KappaName_);

        vectorField n = patch().nf();

        return n & KappaWall & n;
    }
    else
    {
        FatalErrorIn
        (
            "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::Kappa()"
            " const"
        )   << "Did not find field " << KappaName_
            << " on mesh " << mesh.name() << " patch " << patch().name()
            << endl
            << "Please set 'Kappa' to 'none', a valid volScalarField"
            << " or a valid volSymmTensorField." << exit(FatalError);

        return scalarField(0);
    }
}


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the directMappedPatchBase
    const directMappedPatchBase& mpp = refCast<const directMappedPatchBase>
    (
        patch().patch()
    );
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch = refCast<const fvMesh>
    (
        nbrMesh
    ).boundary()[mpp.samplePolyPatch().index()];

    // Force recalculation of mapping and schedule
    const mapDistribute& distMap = mpp.map();

    tmp<scalarField> intFld = patchInternalField();


    // Calculate the temperature by harmonic averaging
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField& nbrField =
    refCast
    <
        const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
    >
    (
        nbrPatch.lookupPatchField<volScalarField, scalar>
        (
            neighbourFieldName_
        )
    );

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld = nbrField.patchInternalField();
    mapDistribute::distribute
    (
        static_cast<Pstream::commsTypes>(Pstream::defaultCommsType()),
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nbrIntFld
    );

    // Swap to obtain full local values of neighbour Kappa*delta
    scalarField nbrKappaDelta = nbrField.Kappa()*nbrPatch.deltaCoeffs();
    mapDistribute::distribute
    (
        static_cast<Pstream::commsTypes>(Pstream::defaultCommsType()),
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nbrKappaDelta
    );

    tmp<scalarField> myKappaDelta = Kappa()*patch().deltaCoeffs();


    // Both sides agree on
    // - temperature : (myKappaDelta*fld + nbrKappaDelta*nbrFld)/(myKappaDelta+nbrKappaDelta)
    // - gradient    : (temperature-fld)*delta
    // We've got a degree of freedom in how to implement this in a mixed bc.
    // (what gradient, what fixedValue and mixing coefficient)
    // Two reasonable choices:
    // 1. specify above temperature on one side (preferentially the high side)
    //    and above gradient on the other. So this will switch between pure
    //    fixedvalue and pure fixedgradient
    // 2. specify gradient and temperature such that the equations are the
    //    same on both sides. This leads to the choice of
    //    - refGradient = zero gradient
    //    - refValue = neighbour value
    //    - mixFraction = nbrKappaDelta / (nbrKappaDelta + myKappaDelta())


    this->refValue() = nbrIntFld;

    this->refGrad() = 0.0;

    this->valueFraction() = nbrKappaDelta / (nbrKappaDelta + myKappaDelta());

    mixedFvPatchScalarField::updateCoeffs();


    if (debug)
    {
        scalar Q = gSum(Kappa()*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " -> "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->dimensionedInternalField().name() << " :"
            << " heatFlux:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("neighbourFieldName")<< neighbourFieldName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("Kappa") << KappaName_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
