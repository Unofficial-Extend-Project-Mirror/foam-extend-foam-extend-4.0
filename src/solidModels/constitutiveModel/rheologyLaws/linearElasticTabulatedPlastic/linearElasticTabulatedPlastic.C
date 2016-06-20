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
    (rheologyLaw, linearElasticTabulatedPlastic, dictionary);
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
    numDiffDelta_(1e-6),
    sigmaY_
    (
        dimensionedScalar
        ("initialYieldStress", dimPressure, stressPlasticStrainSeries_(0.0))
        )
{}


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


//- Return yield stress - cellID needed for multiMaterial
Foam::scalar Foam::linearElasticTabulatedPlastic::
sigmaY(const scalar epsilonPEq, const label cellID) const
{
  return stressPlasticStrainSeries_(epsilonPEq);
}


Foam::tmp<Foam::volScalarField> Foam::linearElasticTabulatedPlastic::Ep() const
{
  // not used in aravasMises
  notImplemented("linearElasticTabulatedPlastic::Ep() is not implemented");

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
            dimensionedScalar("zeroEp", dimPressure, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::linearElasticTabulatedPlastic::Ep(const volScalarField& epsilonPEq) const
{
  // not used in aravasMises
  notImplemented("linearElasticTabulatedPlastic::Ep() is not implemented");

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
            dimensionedScalar("zeroEp", dimPressure, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& Ep = tresult();

    // The user specifies a table of "Stress vs plasticStrain"
    // However, the plasticity return algorithm - in plasticityModel.correct() -
    // requires the plastic modulus at the current Mises stress
    // Therefore, we need essentially need to differentiate the
    // stressVsPlasticStrain curve. There are many ways we can do this,
    // such as pre-smoothing the data
    // and then fitting least squares polynomials.
    // For now, we will keep it simple and assume that the input data is smooth
    // and that the strain increment is relatively small (at least 10 points
    // per 1% strain).
    // We will then find the plastic modelus (Ep) by making a linear fit
    // to the sigmaEq at + and - numDiffDelta strain,
    // where numDiffDelta will be relatively small, depending on resolution.

    const scalarField& epsilonPEqI = epsilonPEq.internalField();
    scalarField& EpI = Ep.internalField();

    // Correction: the current point (epsilonPEq, sigmaEq) may not lie exactly
    // on the yield surface, so we will apply a correction to make sure it
    // returns exactly to it
    //const volScalarField& sigmaEq =
    //mesh().objectRegistry::lookupObject<volScalarField>("sigmaEq");
    //const volScalarField& oldSigmaEq =
      //mesh().objectRegistry::lookupObject<volScalarField>("oldSigmaEq");
    //const volScalarField& oldEpsilonPEq =
    //sigma().mesh().objectRegistry::lookupObject<volScalarField>
    // ("oldEpsilonPEq");
    //const scalarField& sigmaEqI = sigmaEq.internalField();
    //const scalarField& oldSigmaEqI = oldSigmaEq.internalField();
    //const scalarField& oldEpsilonPEqI = oldEpsilonPEq.internalField();


    forAll(epsilonPEqI, celli)
    {
      const scalar epsilonPEqi = epsilonPEq.internalField()[celli];
      // current yield stress for given plastic strain
      //const scalar sigmaY = stressPlasticStrainSeries_(epsilonPEqi);

      // find the yield stress at a slightly higher and lower strain and use
      // these value to calculate plastic modulus
      // This is a very basic method to calculate the modulus but is OK for now
      // assuming the data is relatively smooth
      const scalar lowerEpsilonPEq = max(0.0,epsilonPEqi-numDiffDelta_);
      const scalar higherEpsilonPEq = epsilonPEqi+numDiffDelta_;
      const scalar lowerSigmaY = stressPlasticStrainSeries_(lowerEpsilonPEq);
      const scalar higherSigmaY = stressPlasticStrainSeries_(higherEpsilonPEq);

      // use actual sigmaEq value to correct any deviation from the
      // yield surface
      // const scalar lowerEpsilonPEq = oldEpsilonPEqI[celli];
      // const scalar higherEpsilonPEq = epsilonPEqi;
      // const scalar lowerSigmaY = oldSigmaEqI[celli];
      // const scalar higherSigmaY =
      //stressPlasticStrainSeries_(higherEpsilonPEq);

      // Set plastic modulus
      EpI[celli] =
          (higherSigmaY - lowerSigmaY)
          /max(SMALL, higherEpsilonPEq-lowerEpsilonPEq);
    }

    Ep.correctBoundaryConditions();

    return tresult;
}

// ************************************************************************* //
