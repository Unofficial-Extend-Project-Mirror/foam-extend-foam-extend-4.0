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
    Zeljko Tukovic. FSB Zagreb

-------------------------------------------------------------------------------
*/

#include "expElasticPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(expElasticPlastic, 0);
    addToRunTimeSelectionTable(rheologyLaw, expElasticPlastic, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::expElasticPlastic::expElasticPlastic
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
    sigmaYInf_(dict.lookup("sigmaYInf")),
    delta_(dict.lookup("delta"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::expElasticPlastic::~expElasticPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::expElasticPlastic::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rhoLaw",
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


Foam::tmp<Foam::volScalarField> Foam::expElasticPlastic::E() const
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


Foam::tmp<Foam::volScalarField> Foam::expElasticPlastic::nu() const
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


Foam::tmp<Foam::volScalarField> Foam::expElasticPlastic::sigmaY() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "lawSigmaY",
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


Foam::scalar Foam::expElasticPlastic::sigmaY
(
    const scalar epsilonPEq,
    const label cellID
) const
{
    scalar ePEq = epsilonPEq;

    if (ePEq < 0)
    {
        Info << "Limiting epsilonPEq in expElasticPlastic::sigmaY(...)" << endl;
        ePEq = mag(ePEq);
    }

    // return yield stress given epsilonPEq
    return sigmaY_.value() + ePEq*Ep_.value()
      + (sigmaYInf_.value() - sigmaY_.value())
       *(1.0 - exp(-delta_.value()*ePEq));
}


Foam::scalar Foam::expElasticPlastic::dSigmaY
(
    const scalar epsilonPEq,
    const label cellID
) const
{
    scalar ePEq = epsilonPEq;

    if (ePEq < 0)
    {
        Info << "Limiting epsilonPEq in expElasticPlastic::dSigmaY(...)"
            << endl;
        ePEq = mag(ePEq);
    }

    // return yield stress given epsilonPEq
    return Ep_.value()
      + (sigmaYInf_.value() - sigmaY_.value())
       *delta_.value()*exp(-delta_.value()*ePEq);
}


Foam::tmp<Foam::volScalarField> Foam::expElasticPlastic::Ep() const
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

// ************************************************************************* //
