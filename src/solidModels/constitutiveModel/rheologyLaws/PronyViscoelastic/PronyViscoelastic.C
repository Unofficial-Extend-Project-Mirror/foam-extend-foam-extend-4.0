/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "PronyViscoelastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PronyViscoelastic, 0);
    addToRunTimeSelectionTable(rheologyLaw, PronyViscoelastic, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::PronyViscoelastic::PronyViscoelastic
(
    const word& name,
    const volSymmTensorField& sigma,
    const dictionary& dict
)
:
    rheologyLaw(name, sigma, dict),
    rho_(dict.lookup("rho")),
    k_("k", dict, readInt(dict.lookup("size"))),
    kDimensions_(dict.lookup("kDimensions")),
    tau_("tau", dict, readInt(dict.lookup("size"))),
    tauDimensions_(dict.lookup("tauDimensions")),
    nu_(dict.lookup("nu"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PronyViscoelastic::~PronyViscoelastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::PronyViscoelastic::rho(scalar t) const
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


Foam::tmp<Foam::volScalarField> Foam::PronyViscoelastic::E(scalar t) const
{
    scalar E = 0.0;

    E = k_[0];

    for (int i=1; i<k_.size(); i++)
    {
        E += k_[i]*exp(-t/tau_[i]);
    }

    if (t < 0)
    {
        E = 0;
    }

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
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("E", kDimensions_, E),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();
    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::PronyViscoelastic::nu(scalar t) const
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


Foam::tmp<Foam::volScalarField> Foam::PronyViscoelastic::J(scalar t) const
{
    notImplemented(type() + "::J(scalar t)");

    return 1.0/E(t);
}


Foam::tmp<Foam::volScalarField> Foam::PronyViscoelastic::Ep() const
{
    notImplemented(type() + "::Ep()");

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
            dimensionedScalar("zeroEp", dimForce/dimArea, GREAT),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::PronyViscoelastic::sigmaY() const
{
    notImplemented(type() + "::sigmaY()");

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
            dimensionedScalar("zeroSigmaY", dimForce/dimArea, GREAT),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

// ************************************************************************* //
