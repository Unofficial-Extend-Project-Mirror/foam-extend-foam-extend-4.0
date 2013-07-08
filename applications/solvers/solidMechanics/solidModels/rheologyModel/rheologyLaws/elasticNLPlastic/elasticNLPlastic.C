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

#include "elasticNLPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(elasticNLPlastic, 0);
    addToRunTimeSelectionTable(rheologyLaw, elasticNLPlastic, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::elasticNLPlastic::elasticNLPlastic
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

Foam::elasticNLPlastic::~elasticNLPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::elasticNLPlastic::rho() const
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


Foam::tmp<Foam::volScalarField> Foam::elasticNLPlastic::E() const
{
    // Correction of modulus of elasticity to account for 
    // stress-strain curve continuity!
    // Done according to yield stress value - E = sigmaY/epsY! 

    dimensionedScalar Ecorr = 
        sigmaY_
       /::exp
        (
            log
            (
                log(matStrength_/(matStrength_ - sigmaY_)).value()
               /bCf_.value()
            )
           /nCf_.value()
        ); 

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
            Ecorr,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::elasticNLPlastic::nu() const
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


Foam::tmp<Foam::volScalarField> Foam::elasticNLPlastic::sigmaY() const
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


Foam::tmp<Foam::volScalarField> Foam::elasticNLPlastic::Ep() const
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


Foam::tmp<Foam::volScalarField> Foam::elasticNLPlastic::
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

    const scalarField& sigmaEqI = sigmaEq.internalField();

    scalar epsY =  exp ( log ( log(matStrength_/(matStrength_ -
        sigmaY_)).value() /bCf_.value() ) /nCf_.value() );

    dimensionedScalar Ecorr =  sigmaY_ /epsY; 

    forAll(sigmaEqI, cellI)
    {

      scalar epsCurrI =  exp ( log ( log(matStrength_.value()/(max(matStrength_.value()/1e6,matStrength_.value() - sigmaEqI[cellI]))) /bCf_.value()) /nCf_.value());

      dimensionedScalar Ep = matStrength_*bCf_*nCf_ *pow(epsCurrI, nCf_ - 1.0)
        *exp(-bCf_*pow(epsCurrI, nCf_));

      tresult().internalField()[cellI] = 
            Ep.value()/(1.0 - Ep.value()/Ecorr.value());
    }

    tresult().correctBoundaryConditions();

    return tresult;
}


// ************************************************************************* //
