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

#include "turbulentTemperatureCoupledBaffleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "directMappedPatchBase.H"
#include "regionProperties.H"
#include "basicThermo.H"
#include "RASModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::interfaceOwner
(
    const polyMesh& nbrRegion,
    const polyPatch& nbrPatch
) const
{
    const fvMesh& myRegion = patch().boundaryMesh().mesh();

    if (nbrRegion.name() == myRegion.name())
    {
        return patch().index() < nbrPatch.index();
    }
    else
    {
        const regionProperties& props =
            myRegion.objectRegistry::parent().lookupObject<regionProperties>
            (
                "regionProperties"
            );

        label myIndex = findIndex(props.fluidRegionNames(), myRegion.name());
        if (myIndex == -1)
        {
            label i = findIndex(props.solidRegionNames(), myRegion.name());

            if (i == -1)
            {
                FatalErrorIn
                (
                    "turbulentTemperatureCoupledBaffleFvPatchScalarField"
                    "::interfaceOwner(const polyMesh&"
                    ", const polyPatch&)const"
                )   << "Cannot find region " << myRegion.name()
                    << " neither in fluids " << props.fluidRegionNames()
                    << " nor in solids " << props.solidRegionNames()
                    << exit(FatalError);
            }
            myIndex = props.fluidRegionNames().size() + i;
        }
        label nbrIndex = findIndex
        (
            props.fluidRegionNames(),
            nbrRegion.name()
        );
        if (nbrIndex == -1)
        {
            label i = findIndex(props.solidRegionNames(), nbrRegion.name());

            if (i == -1)
            {
                FatalErrorIn
                (
                    "coupleManager::interfaceOwner"
                    "(const polyMesh&, const polyPatch&) const"
                )   << "Cannot find region " << nbrRegion.name()
                    << " neither in fluids " << props.fluidRegionNames()
                    << " nor in solids " << props.solidRegionNames()
                    << exit(FatalError);
            }
            nbrIndex = props.fluidRegionNames().size() + i;
        }

        return myIndex < nbrIndex;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::
turbulentTemperatureCoupledBaffleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    neighbourFieldName_("undefined-neighbourFieldName"),
    KappaName_("undefined-Kappa")
{}


Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::
turbulentTemperatureCoupledBaffleFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    neighbourFieldName_(ptf.neighbourFieldName_),
    KappaName_(ptf.KappaName_)
{}


Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::
turbulentTemperatureCoupledBaffleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    neighbourFieldName_(dict.lookup("neighbourFieldName")),
    KappaName_(dict.lookup("Kappa"))
{
    if (!isA<directMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "turbulentTemperatureCoupledBaffleFvPatchScalarField::"
            "turbulentTemperatureCoupledBaffleFvPatchScalarField\n"
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
}


Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::
turbulentTemperatureCoupledBaffleFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wtcsf, iF),
    neighbourFieldName_(wtcsf.neighbourFieldName_),
    KappaName_(wtcsf.KappaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::Kappa() const
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
            "turbulentTemperatureCoupledBaffleFvPatchScalarField::Kappa() const"
        )   << "Did not find field " << KappaName_
            << " on mesh " << mesh.name() << " patch " << patch().name()
            << endl
            << "Please set 'Kappa' to 'none', a valid volScalarField"
            << " or a valid volSymmTensorField." << exit(FatalError);

        return scalarField(0);
    }
}


void Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::updateCoeffs()
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
    (void)distMap.schedule();

    tmp<scalarField> intFld = patchInternalField();

    if (interfaceOwner(nbrMesh, nbrPatch.patch()))
    {
        // Note: other side information could be cached - it only needs
        // to be updated the first time round the iteration (i.e. when
        // switching regions) but unfortunately we don't have this information.


        // Calculate the temperature by harmonic averaging
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        const turbulentTemperatureCoupledBaffleFvPatchScalarField& nbrField =
        refCast<const turbulentTemperatureCoupledBaffleFvPatchScalarField>
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

        // Calculate common wall temperature. Reuse *this to store common value.
        scalarField Twall
        (
            (myKappaDelta()*intFld() + nbrKappaDelta*nbrIntFld)
          / (myKappaDelta() + nbrKappaDelta)
        );
        // Assign to me
        fvPatchScalarField::operator=(Twall);
        // Distribute back and assign to neighbour
        mapDistribute::distribute
        (
            static_cast<Pstream::commsTypes>(Pstream::defaultCommsType()),
            distMap.schedule(),
            nbrField.size(),
            distMap.constructMap(),     // reverse : what to send
            distMap.subMap(),
            Twall
        );
        const_cast<turbulentTemperatureCoupledBaffleFvPatchScalarField&>
        (
            nbrField
        ).fvPatchScalarField::operator=(Twall);
    }

    if (debug)
    {
        //tmp<scalarField> normalGradient =
        //    (*this-intFld())
        //  * patch().deltaCoeffs();

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

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedValueFvPatchScalarField::write(os);
    os.writeKeyword("neighbourFieldName")<< neighbourFieldName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("Kappa") << KappaName_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureCoupledBaffleFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
