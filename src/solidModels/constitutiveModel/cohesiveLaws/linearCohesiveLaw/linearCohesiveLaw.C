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

\*---------------------------------------------------------------------------*/

#include "linearCohesiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"
#include "cohesiveFvPatch.H"
#include "cohesiveLaw.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearCohesiveLaw, 0);
    addToRunTimeSelectionTable(cohesiveLaw, linearCohesiveLaw, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearCohesiveLaw::linearCohesiveLaw
(
    const word& name,
    const volSymmTensorField& sigma,
    const dictionary& dict
)
:
  cohesiveLaw(name, sigma, dict),
  GIc_(dict.lookup("GIc")),
  GIIc_(dict.lookup("GIIc")),
  sigmaMax_(dict.lookup("sigmaMax")),
  tauMax_(dict.lookup("tauMax"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearCohesiveLaw::~linearCohesiveLaw()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::linearCohesiveLaw::materials() const
{
  notImplemented(type() + "::materials()");

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "materials",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
}

Foam::tmp<Foam::surfaceScalarField> Foam::linearCohesiveLaw::sigmaMax() const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "sigmaMax",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            sigmaMax_
        )
    );
}

Foam::tmp<Foam::surfaceScalarField> Foam::linearCohesiveLaw::tauMax() const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "tauMax",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            tauMax_
        )
    );
}

Foam::tmp<Foam::surfaceScalarField> Foam::linearCohesiveLaw::GIc() const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "GIc",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            GIc_
        )
    );
}

Foam::tmp<Foam::surfaceScalarField> Foam::linearCohesiveLaw::GIIc() const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "GIIc",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            GIIc_
        )
    );
}


void Foam::linearCohesiveLaw::damageTractions
(
 scalar& tN,
 scalar& tS,
 const scalar deltaN,
 const scalar deltaS,
 const scalar GI,
 const scalar GII,
 const label faceID,
 const scalarField& globalPatchMaterials
 ) const
{
  // Step 1
  // Calculate apparent initiation normal traction tNia
  // tNia is the initiation traction assuming that the ratio tN/tS
  // has stayed constant at the current value i.e.
  // tNia = B * tN
  // tSia = B * tS
  // where we term B the initiation traction constant.

  // Calculate current B
  const scalar B = 1.0 /
    ::sqrt(
       (tN/sigmaMax_.value())*(tN/sigmaMax_.value())
       + (tS/tauMax_.value())*(tS/tauMax_.value())
       );

  // Calculate apparent initiation tractions
  const scalar tNia = B * tN;
  const scalar tSia = B * tS;

  // Step 2
  // Calculate apparent critical delta deltaCa
  // deltaCa would be the deltaC if the mode mix stayed constant at
  // the current value
  // GIca = C * GII
  // GIIca = C * GII
  // where we term C the mode-mix constant

  // Calculate C
  const scalar C = 1.0 /
    ( (GI/GIc_.value()) + (GII/GIIc_.value()) );

  // Calculate apparent critical energies
  const scalar GIca = C * GI;
  const scalar GIIca = C * GII;

  // Calculate apparent critical deltas
  const scalar deltaNca = 2.0 * GIca / tNia;
  const scalar deltaSca = 2.0 * GIIca / tSia;

  // Step 3
  // Set current tractions linearly decreasing from
  // tia to deltaCa
  // When the face is close to propagation the tN and
  // tS will get close to zero so we will limit
  // them to be positive as they are magnitudes
  tN = max( 0.0, tNia * ( 1.0 - (deltaN/deltaNca) ) );
  tS = max( 0.0, tSia * ( 1.0 - (deltaS/deltaSca) ) );

  // Step 4
  // Hold the effective traction traction and
  // allow the mix of tractions to vary based on:
  // (tN/tS) = a * (deltaN/deltaS)
  // where a=1 is assumed here.
  // This essentially allows the mode-mixity
  // to vary depending on the deltas
  tN = (
    (sigmaMax_.value()/B) * deltaN /
    (SMALL + ::sqrt(
            (deltaN*deltaN)
            + (deltaS*deltaS)*(sigmaMax_.value()*sigmaMax_.value()
                               /(tauMax_.value()*tauMax_.value()))
            ))
    );
  tS = (
    (tauMax_.value()/B) * deltaS /
    (SMALL + ::sqrt(
            (deltaS*deltaS)
            + (deltaN*deltaN)*(tauMax_.value()*tauMax_.value()
                               /(sigmaMax_.value()*sigmaMax_.value()))
            ))
    );
}


Foam::tmp<Foam::surfaceVectorField>
Foam::linearCohesiveLaw::interfaceTraction
(
 surfaceVectorField n,
 volVectorField U,
 volTensorField gradU,
 volScalarField mu,
 volScalarField lambda
 ) const
{
  notImplemented(type() + "::interfaceTraction()");

    tmp<surfaceVectorField> tresult
    (
        new surfaceVectorField
        (
            IOobject
            (
            "interfaceTraction",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
            mesh(),
            dimensionedVector("zero", dimForce/dimArea, vector(0, 0, 0))
    )
    );

    return tresult;
}

// ************************************************************************* //
