// The FOAM Project // File: elasticPlastic.C
/*
-------------------------------------------------------------------------------
 =========         | Class Implementation
 \\      /         |
  \\    /          | Name:   elasticPlastic
   \\  /           | Family: rheologyLaw
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2004 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
DESCRIPTION

AUTHOR
    Hrvoje Jasak.

-------------------------------------------------------------------------------
*/

#include "nonLinearElasticPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonLinearElasticPlastic, 0);
    addToRunTimeSelectionTable(rheologyLaw, nonLinearElasticPlastic, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::nonLinearElasticPlastic::nonLinearElasticPlastic
(
    const word& name,
    const volSymmTensorField& sigma,
    const dictionary& dict
)
:
    rheologyLaw(name, sigma, dict),
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    sigmaY_(dict.lookup("sigmaY")),
    Ep_(dict.lookup("Ep")),
    matStrength_(dict.lookup("sigmaMax")),
    bCf_(dict.lookup("bCf")),
    nCf_(dict.lookup("nCf"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonLinearElasticPlastic::~nonLinearElasticPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::nonLinearElasticPlastic::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::nonLinearElasticPlastic::E() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            E_,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::nonLinearElasticPlastic::
E(const volScalarField& epsEq) const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zeroE", dimPressure, E_.value()),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    
    scalar epsY =  exp ( log ( log(matStrength_/(matStrength_ -
    sigmaY_)).value() /bCf_.value() ) /nCf_.value() );

    const scalarField& epsEqI = epsEq.internalField();

    forAll(epsEqI, cellI)
    {

        dimensionedScalar E =  matStrength_*bCf_*nCf_ *pow(epsEqI[cellI], nCf_
            - 1.0) *exp(-bCf_*pow(epsEqI[cellI], nCf_));

        // Correction of initial modulus to avoid infinity/GREAT 
        // for small strains and for unloading
	// strain of 0.1% might be wrong for some materials

        if(epsEqI[cellI] < 0.001)
        {
            E = matStrength_*bCf_*nCf_*pow(0.001, nCf_ - 1.0)
               *exp(-bCf_*pow(0.001, nCf_));
        }

        if(epsEqI[cellI] > epsY)
        { 
            E = matStrength_*bCf_*nCf_*pow(epsY, nCf_ - 1.0)
               *exp(-bCf_*pow(epsY, nCf_));
        }

        tresult().internalField()[cellI] = E.value();
    }

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::nonLinearElasticPlastic::nu() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "nu",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            nu_,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::nonLinearElasticPlastic::sigmaY() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "sigmaY",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            sigmaY_,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::nonLinearElasticPlastic::Ep() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Ep",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            Ep_,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::nonLinearElasticPlastic::
Ep(const volScalarField& sigmaEq) const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "Ep",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zeroEp", dimPressure, Ep_.value()),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    scalar epsY =  exp ( log ( log(matStrength_/(matStrength_ -
    sigmaY_)).value() /bCf_.value() ) /nCf_.value() );

    dimensionedScalar E =  matStrength_*bCf_*nCf_ *pow(epsY, nCf_ - 1.0)
        *exp(-bCf_*pow(epsY, nCf_));

    const scalarField& sigmaEqI = sigmaEq.internalField();

    forAll(sigmaEqI, cellI)
    {

      scalar epsCurrI =  exp ( log ( log(matStrength_.value()/(max(matStrength_.value()/1e6,matStrength_.value() - sigmaEqI[cellI]))) /bCf_.value()) /nCf_.value());

      dimensionedScalar Ep = matStrength_*bCf_*nCf_ *pow(epsCurrI, nCf_ - 1.0)
        *exp(-bCf_*pow(epsCurrI, nCf_));

        tresult().internalField()[cellI] = 
            Ep.value()/(1.0 - Ep.value()/E.value());
    }

    tresult().correctBoundaryConditions();

    return tresult;
}


// ************************************************************************* //
