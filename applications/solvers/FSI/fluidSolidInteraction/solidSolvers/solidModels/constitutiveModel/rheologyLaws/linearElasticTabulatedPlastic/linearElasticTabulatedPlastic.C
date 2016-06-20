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

Class
    linearElasticTabulatedPlastic

\*---------------------------------------------------------------------------*/

#include "linearElasticTabulatedPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElasticTabulatedPlastic, 0);
    addToRunTimeSelectionTable
    (
        rheologyLaw, linearElasticTabulatedPlastic, dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearElasticTabulatedPlastic::linearElasticTabulatedPlastic
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
    stressPlasticStrainSeries_(dict),
    sigmaY_
    (
        dimensionedScalar
        (
            "initialYieldStress", dimPressure, stressPlasticStrainSeries_(0.0)
        )
    )
{
    Info << "sigmaY: " << sigmaY_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearElasticTabulatedPlastic::~linearElasticTabulatedPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::linearElasticTabulatedPlastic::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "lawRho",
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


Foam::tmp<Foam::volScalarField> Foam::linearElasticTabulatedPlastic::E() const
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


Foam::tmp<Foam::volScalarField> Foam::linearElasticTabulatedPlastic::nu() const
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


Foam::tmp<Foam::volScalarField>
Foam::linearElasticTabulatedPlastic::sigmaY() const
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
                // IOobject::NO_WRITE
                IOobject::AUTO_WRITE
            ),
            mesh(),
            sigmaY_,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


//- Return yield stress - cellID needed for multiMaterial
Foam::scalar Foam::linearElasticTabulatedPlastic::sigmaY
(
    const scalar epsilonPEq, const label cellID
) const
{
  return stressPlasticStrainSeries_(epsilonPEq);
}


//- Return yield stress derivative - cellID needed for multiMaterial
Foam::scalar Foam::linearElasticTabulatedPlastic::dSigmaY
(
    const scalar epsilonPEq, const label cellID
) const
{
    notImplemented
    (
        "linearElasticTabulatedPlastic::dSigmaY() is not implemented"
    );

    return 0;
}


Foam::tmp<Foam::volScalarField> Foam::linearElasticTabulatedPlastic::Ep() const
{
    if (stressPlasticStrainSeries_.size() > 2)
    {
        notImplemented
        (
            "linearElasticTabulatedPlastic::Ep() is not implemented for"
            " nonLinear plasticity"
        );
    }

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
            dimensionedScalar("zeroEp", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (stressPlasticStrainSeries_.size() > 1)
    {
        // Calculate plastic modulus from two points
        tresult() =
            (
                stressPlasticStrainSeries_[1].second()
                - stressPlasticStrainSeries_[0].second()
            )
            /(
                stressPlasticStrainSeries_[1].first()
                - stressPlasticStrainSeries_[0].first()
             );

        tresult().correctBoundaryConditions();
    }

    return tresult;
}


Foam::tmp<Foam::volScalarField>
Foam::linearElasticTabulatedPlastic::Ep(const volScalarField& epsilonPEq) const
{
    notImplemented
    (
        "linearElasticTabulatedPlastic::Ep(const volScalarField& epsilonPEq) "
        "is not implemented"
    );

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
            dimensionedScalar("zeroEp", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


bool Foam::linearElasticTabulatedPlastic::nonLinearPlasticity() const
{
    // if there are more than two points in the yield stress series then
    // the plastic law is nonlinear
    if (stressPlasticStrainSeries_.size() > 2)
    {
        return true;
    }

    return false;
}

// ************************************************************************* //
