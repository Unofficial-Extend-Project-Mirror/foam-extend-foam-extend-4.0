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

Description

\*---------------------------------------------------------------------------*/

#include "MaxwellViscoelastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MaxwellViscoelastic, 0);
    addToRunTimeSelectionTable(rheologyLaw, MaxwellViscoelastic, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::MaxwellViscoelastic::MaxwellViscoelastic
(
    const word& name,
    const volSymmTensorField& sigma,
    const dictionary& dict
)
:
    rheologyLaw(name, sigma, dict),
    rho_(dict.lookup("rho")),
    k_(dict.lookup("k")),
    eta_(dict.lookup("eta")),
    nu_(dict.lookup("nu"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MaxwellViscoelastic::~MaxwellViscoelastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::MaxwellViscoelastic::rho(scalar t) const
{
    tmp<volScalarField> tresult
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

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::MaxwellViscoelastic::E(scalar t) const
{
    scalar tau = eta_.value()/k_.value();

    tmp<volScalarField> tE
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
            k_*exp(-t/tau),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (t < 0)
    {
        tE().internalField() = 0.0;
        tE().correctBoundaryConditions();
    }

    return tE;
}


Foam::tmp<Foam::volScalarField> Foam::MaxwellViscoelastic::nu(scalar t) const
{
    tmp<volScalarField> tresult
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

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::MaxwellViscoelastic::J(scalar t) const
{
    tmp<volScalarField> tJ
    (
        new volScalarField
        (
            IOobject
            (
                "J",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar
            (
                "J",
                dimless/k_.dimensions(),
                1.0/k_.value() + t/eta_.value()
            ),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (t < 0)
    {
        tJ().internalField() = 0.0;
        tJ().correctBoundaryConditions();
    }

    return tJ;
}


// ************************************************************************* //
