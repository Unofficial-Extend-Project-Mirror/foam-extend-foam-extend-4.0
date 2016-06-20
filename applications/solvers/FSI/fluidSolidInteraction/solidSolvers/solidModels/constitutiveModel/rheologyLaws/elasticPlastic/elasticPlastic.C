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

#include "elasticPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(elasticPlastic, 0);
    addToRunTimeSelectionTable(rheologyLaw, elasticPlastic, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::elasticPlastic::elasticPlastic
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
    Ep_(dict.lookup("Ep"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::elasticPlastic::~elasticPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::elasticPlastic::rho() const
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


Foam::tmp<Foam::volScalarField> Foam::elasticPlastic::E() const
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


Foam::tmp<Foam::volScalarField> Foam::elasticPlastic::nu() const
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


Foam::tmp<Foam::volScalarField> Foam::elasticPlastic::sigmaY() const
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
                IOobject::AUTO_WRITE
            ),
            mesh(),
            sigmaY_,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

Foam::scalar Foam::elasticPlastic::sigmaY(const scalar epsilonPEq, const label cellID) const
{
  // return yield stress given epsilonPEq
  return sigmaY_.value() + epsilonPEq*Ep_.value();
}


Foam::tmp<Foam::volScalarField> Foam::elasticPlastic::Ep() const
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
