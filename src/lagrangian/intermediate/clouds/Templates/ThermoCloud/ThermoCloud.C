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

#include "ThermoCloud.H"
#include "HeatTransferModel.H"

#include "interpolationCellPoint.H"
#include "ThermoParcel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ThermoCloud<ParcelType>::ThermoCloud
(
    const word& cloudType,
    const volPointInterpolation& vpi,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    basicThermo& thermo
)
:
    KinematicCloud<ParcelType>
    (
        cloudType,
        vpi,
        rho,
        U,
        thermo.mu(),
        g
    ),
    thermoCloud(),
    constProps_(this->particleProperties()),
    carrierThermo_(thermo),
    heatTransferModel_
    (
        HeatTransferModel<ThermoCloud<ParcelType> >::New
        (
            this->particleProperties(),
            *this
        )
    ),
    TIntegrator_
    (
        scalarIntegrationScheme::New
        (
            "T",
            this->particleProperties().subDict("integrationSchemes")
        )
    ),
    radiation_(this->particleProperties().lookup("radiation")),
    hTrans_
    (
        IOobject
        (
            this->name() + "hTrans",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh(),
        dimensionedScalar("zero", dimensionSet(1, 2, -2, 0, 0), 0.0)
    ),
    hCoeff_
    (
        IOobject
        (
            this->name() + "hCoeff",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh(),
        dimensionedScalar("zero", dimensionSet(1, 2, -3, -1, 0), 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ThermoCloud<ParcelType>::~ThermoCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::ThermoCloud<ParcelType>::resetSourceTerms()
{
    KinematicCloud<ParcelType>::resetSourceTerms();
    hTrans_.field() = 0.0;
    hCoeff_.field() = 0.0;
}


template<class ParcelType>
void Foam::ThermoCloud<ParcelType>::evolve()
{
    const volScalarField& T = carrierThermo_.T();
    const volScalarField cp = carrierThermo_.Cp();

    autoPtr<interpolation<scalar> > rhoInterpolator = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->vpi(),
        this->rho()
    );

    autoPtr<interpolation<vector> > UInterpolator = interpolation<vector>::New
    (
        this->interpolationSchemes(),
        this->vpi(),
        this->U()
    );

    autoPtr<interpolation<scalar> > muInterpolator = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->vpi(),
        this->mu()
    );

    autoPtr<interpolation<scalar> > TInterpolator = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->vpi(),
        T
    );

    autoPtr<interpolation<scalar> > cpInterpolator = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->vpi(),
        cp
    );

    typename ParcelType::trackData td
    (
        *this,
        constProps_,
        rhoInterpolator(),
        UInterpolator(),
        muInterpolator(),
        TInterpolator(),
        cpInterpolator(),
        this->g().value()
    );

    inject(td);

    this->move(td);
}


template<class ParcelType>
template<class TrackingData>
void Foam::ThermoCloud<ParcelType>::inject
(
    TrackingData& td
)
{
    // Injection is same as for KinematicCloud<ParcelType>
    KinematicCloud<ParcelType>::inject(td);
}


// ************************************************************************* //
