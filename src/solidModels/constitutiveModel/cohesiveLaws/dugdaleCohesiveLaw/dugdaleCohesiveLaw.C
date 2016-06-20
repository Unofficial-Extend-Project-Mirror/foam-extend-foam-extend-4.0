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

#include "dugdaleCohesiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"
#include "cohesiveFvPatch.H"
#include "cohesiveLaw.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dugdaleCohesiveLaw, 0);
    addToRunTimeSelectionTable(cohesiveLaw, dugdaleCohesiveLaw, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::dugdaleCohesiveLaw::dugdaleCohesiveLaw
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

Foam::dugdaleCohesiveLaw::~dugdaleCohesiveLaw()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dugdaleCohesiveLaw::materials() const
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

Foam::tmp<Foam::surfaceScalarField> Foam::dugdaleCohesiveLaw::sigmaMax() const
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

Foam::tmp<Foam::surfaceScalarField> Foam::dugdaleCohesiveLaw::tauMax() const
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

Foam::tmp<Foam::surfaceScalarField> Foam::dugdaleCohesiveLaw::GIc() const
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

Foam::tmp<Foam::surfaceScalarField> Foam::dugdaleCohesiveLaw::GIIc() const
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


void Foam::dugdaleCohesiveLaw::damageTractions
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
  // Update normal traction
  tN = (
    sigmaMax_.value() * deltaN /
    (SMALL + ::sqrt( (deltaN*deltaN)
                     + (deltaS*deltaS)*(sigmaMax_.value()*sigmaMax_.value()
                                        /(tauMax_.value()*tauMax_.value()))
        )
        )
    );

  // Update shear traction
  tS = (
    tauMax_.value() * deltaS /
    (SMALL + ::sqrt( (deltaS*deltaS)
                     + (deltaN*deltaN)*(tauMax_.value()*tauMax_.value()
                                        /(sigmaMax_.value()*sigmaMax_.value()))
        )
        )
    );
}


Foam::tmp<Foam::surfaceVectorField>
Foam::dugdaleCohesiveLaw::interfaceTraction
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
