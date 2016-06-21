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

#include "aravasMises.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "constitutiveModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(aravasMises, 0);
    addToRunTimeSelectionTable(plasticityStressReturn, aravasMises, dictionary);


  // static variables
  scalar aravasMises::LoopTol_ = 1e-8;
  label aravasMises::MaxNewtonIter_ = 1000;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
aravasMises::aravasMises
(
 const word& name,
 constitutiveModel& constitutiveModel
)
:
  plasticityStressReturn(name, constitutiveModel),
  constitutiveModel_(constitutiveModel),
  sigmaY_
  (
   IOobject
   (
    "sigmaY",
    constitutiveModel.sigma().time().timeName(),
    constitutiveModel.sigma().db(),
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
    ),
   constitutiveModel.sigmaY()
   ),
  DSigmaY_
  (
   IOobject
   (
    "DSigmaY",
    sigmaY_.time().timeName(),
    sigmaY_.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
   sigmaY_.mesh(),
   dimensionedScalar("0", dimPressure, 0.0)
   ),
  DEpsilonP_
  (
   IOobject
   (
    "DEpsilonP",
    sigmaY_.time().timeName(),
    sigmaY_.db(),
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
    ),
   sigmaY_.mesh(),
   dimensionedSymmTensor("0", dimless, symmTensor::zero)
   ),
  DEpsilonPEq_
  (
   IOobject
   (
    "DEpsilonPEq",
    sigmaY_.time().timeName(),
    sigmaY_.db(),
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
    ),
   sigmaY_.mesh(),
   dimensionedScalar("0", dimless, 0.0)
   ),
  activeYield_
  (
   IOobject
   (
    "activeYield",
    sigmaY_.time().timeName(),
    sigmaY_.db(),
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
    ),
   sigmaY_.mesh(),
   dimensionedScalar("0", dimless, 0)
   ),
  plasticN_
  (
   IOobject
   (
    "plasticN",
    sigmaY_.time().timeName(),
    sigmaY_.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
   sigmaY_.mesh(),
   dimensionedSymmTensor("zero", dimless, symmTensor::zero)
   )
{
  Info << "Creating AravasMises stress return method" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

aravasMises::~aravasMises()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void aravasMises::correct()
{
  // Fields for integration of plasticity model

  // mesh
  const fvMesh& mesh = sigmaY_.mesh();

  // elastic properties
  const volScalarField& mu =
    mesh.objectRegistry::lookupObject<volScalarField>("mu");
  const volScalarField& lambda =
    mesh.objectRegistry::lookupObject<volScalarField>("lambda");

  // Old plastic strain
  const volSymmTensorField& oldEpsilonP =
    mesh.objectRegistry::lookupObject<volSymmTensorField>("epsilonP");

  // Old plastic equivalent strain
  const volScalarField oldEpsilonPEq = sqrt((2.0/3.0)*magSqr(dev(oldEpsilonP)));

  // Current strain increment - should be recalculated within momentum loop
  // before rheology.correct()
  const volSymmTensorField& DEpsilon =
    mesh.objectRegistry::lookupObject<volSymmTensorField>("DEpsilon");

  // Old total strain
  const volSymmTensorField& oldEpsilon =
    mesh.objectRegistry::lookupObject<volSymmTensorField>("epsilon");

  // New total strain
  const volSymmTensorField newEpsilon = oldEpsilon + DEpsilon;

  // Magnitude of current strain increment
  volScalarField magDEpsilon = Foam::magSqr(DEpsilon);

  // Calculate Old Elastic Strain
  volSymmTensorField oldEpsilonElastic = oldEpsilon - oldEpsilonP;

  // Calculate old Elastic Stress
  volSymmTensorField oldSigmaElastic =
      2.0*mu*oldEpsilonElastic + lambda*tr(oldEpsilonElastic)*I;

  // Calculate old Elastic Stress deviatoric component
  volSymmTensorField oldSigmaDev = dev(oldSigmaElastic);

  // Calculate old equivalent stress
  volScalarField oldSigmaEqElastic = sqrt((3.0/2.0)*Foam::magSqr(oldSigmaDev));

  // Calculate new Elastic Strain - this is a trial elastic stress
  volSymmTensorField newEpsilonElastic = newEpsilon - oldEpsilonP;

  // Calculate new Elastic Stress
  volSymmTensorField newSigmaElastic =
      2.0*mu*newEpsilonElastic + lambda*tr(newEpsilonElastic)*I;

  // Calculate change in elastic Stress
  volSymmTensorField DSigmaElastic = newSigmaElastic - oldSigmaElastic;

  // Calculate new elastic Deviatoric
  volSymmTensorField newSigmaDevElastic = dev(newSigmaElastic);

  // Calculate new elastic equivalent stress
  volScalarField newSigmaEqElastic = sqrt((3.0/2.0)*magSqr(newSigmaDevElastic));

  // Store previous iteration for under-relaxation
  DEpsilonP_.storePrevIter();

  // Normalise residual in Newton method with respect to DEpsilon
  const scalar maxMagDEpsilon = max(gMax(mag(DEpsilon.internalField())), SMALL);

  // Calculate beta field (fully elastic-plastic fraction)
  //# include "aravasMisesUpdateBeta.H"

  // Calculate DEpsilonPEq_ and plasticN_
  forAll (DEpsilonP_, cellI)
    {
      // Calculate return direction plasticN
      if (oldSigmaEqElastic[cellI] < SMALL)
      {
          plasticN_[cellI] = I;
      }
      else
      {
          plasticN_[cellI] =
              (3.0/(2.0*newSigmaEqElastic[cellI]))*newSigmaDevElastic[cellI];
      }

      // J2 yield function
      scalar fy =
          Foam::pow((newSigmaEqElastic[cellI]/sigmaY_[cellI]),2.0) - 1.0;

      // Calculate DEpsilonPEq
      if (fy < SMALL)
    {
      // elastic
      DSigmaY_[cellI] = 0.0;
      DEpsilonPEq_[cellI] = 0.0;
    }
      else
    {
      // yielding

      // total equivalent strain matrix (t is start of time-step)
      const scalar ebart = oldEpsilonPEq[cellI];
      // yield stress at start of time-step
      const scalar s00 = sigmaY_[cellI];
      scalar s0 = 0.0; // updated in loop below
      const scalar qe = newSigmaEqElastic[cellI];
      const scalar G = mu[cellI];

      // This Aravas loop calculates deq using Newtons's method
      // Calculate deq and s0
      aravasNewtonLoop
          (DEpsilonPEq_[cellI], s0, s00, G, ebart, qe, cellI, maxMagDEpsilon);

      // Update equivalent plastic strain and increment of yield stress
      DSigmaY_[cellI] = s0 - s00;
    }
    }

  forAll (plasticN_.boundaryField(), patchI)
    {
      if (!plasticN_.boundaryField()[patchI].coupled())
    {
      const labelList& faceCells = mesh.boundary()[patchI].faceCells();
      forAll(plasticN_.boundaryField()[patchI], faceI)
      {
          // Calculate direction plasticN
          if (oldSigmaEqElastic.boundaryField()[patchI][faceI] < SMALL)
          {
              plasticN_.boundaryField()[patchI][faceI] = I;
          }
          else
          {
              plasticN_.boundaryField()[patchI][faceI] =
                  (3.0/(2.0*newSigmaEqElastic.boundaryField()[patchI][faceI]))
                  *newSigmaDevElastic.boundaryField()[patchI][faceI];
          }

          // J2 yield function
          scalar fy =
        Foam::pow(
              (newSigmaEqElastic.boundaryField()[patchI][faceI]
               /
               sigmaY_.boundaryField()[patchI][faceI])
              ,2.0
              )
        - 1.0;

          // Calculate DEpsilonPEq
          if (fy < SMALL)
        {
          // elasticity
          DSigmaY_.boundaryField()[patchI][faceI] = 0.0;
          DEpsilonPEq_.boundaryField()[patchI][faceI] = 0.0;
        }
          else
        {
          // yielding

          // total equivalent strain matrix (t is start of time-step)
          const scalar ebart = oldEpsilonPEq.boundaryField()[patchI][faceI];
          // yield stress at start of time-step
          const scalar s00 = sigmaY_.boundaryField()[patchI][faceI];
          scalar s0 = 0.0; // updated in loop below
          const scalar qe = newSigmaEqElastic.boundaryField()[patchI][faceI];
          const scalar G = mu.boundaryField()[patchI][faceI];

          // Calculate deq and s0
          aravasNewtonLoop
              (
                  DEpsilonPEq_.boundaryField()[patchI][faceI],
                  s0,
                  s00,
                  G,
                  ebart,
                  qe,
                  faceCells[faceI],
                  maxMagDEpsilon
                  );

          // Update increment of yield stress
          DSigmaY_.boundaryField()[patchI][faceI] = s0 - s00;
        }
      }
    }
    }

  // update beta
  // beta is 1.0 for cells that were plastic in the previous step
  // beta is 0.0 for fully elastic cells
  // 0.0 < beta < 1.0 for cells which are crossing the yield surface
  // and beta is used to partition the step into fully elastic and
  // elastic-plastic
  // not correct for Aravas method -> need to check
  //# include "aravasMisesUpdateBeta.H"

  // Update plastic strain increment
  //DEpsilonP_ = beta_*DEpsilonPEq_*plasticN_;
  DEpsilonP_ = DEpsilonPEq_*plasticN_;
  DEpsilonP_.correctBoundaryConditions();

  // Relax plastic strain increment to help convergence
  DEpsilonP_.relax();

  DSigmaY_.correctBoundaryConditions();
}


scalar
aravasMises::qfun (const scalar qe, const scalar deq, const scalar G) const
{
    return (qe - 3*G*deq);
}


// yield stress
scalar aravasMises::s0fun (const scalar ebar, const label cellID) const
{
    scalar ebarNew = ebar;
    if (ebarNew < SMALL)
    {
        ebarNew = SMALL;
    }

    return constitutiveModel_.sigmaY(ebarNew, cellID);
  }

  scalar aravasMises::ebarfun
  (
   const scalar ebart,
   const scalar q,
   const scalar deq,
   const scalar s0
   ) const
  {
    // Optimise integration for large steps
    return (ebart + (q*deq/s0));
  }

  scalar aravasMises::gfun
  (
   const scalar ebart,
   const scalar deq,
   const scalar qe,
   const scalar G,
   const label cellID
   ) const
  {
    const scalar q = qfun(qe, deq, G);

    scalar s0 = s0fun(ebart, cellID);
    const scalar ebar = ebarfun(ebart,q, deq, s0);
    s0 = s0fun(ebar, cellID) ;

    return (Foam::pow(q/s0, 2) - 1.0 );
  }

  void aravasMises::aravasNewtonLoop
  (
   scalar& deq,
   scalar& s0,
   const scalar s00,
   const scalar G,
   const scalar ebart,
   const scalar qe,
   const label cellID,
   const scalar maxMagDEpsilon
   ) const
  {
    // Loop to determine DEpsilonPEq
    // using Newtion's method
    int i = 0;
    //scalar LoopTol = 1e-12;
    //scalar MaxNewtonIter = 1e3;
    scalar fdeq = gfun(ebart, deq, qe, G, cellID);
    // dx is delta for finite difference differentiation
    scalar dx = 0.25e-6;
    scalar residual = 1.0;
    do
      {
    // Numerically calculate parial derivatives using finite differences
    scalar fxplusq  = gfun(ebart, deq+dx, qe, G, cellID);
    scalar fxminusq  = gfun(ebart, deq-dx, qe, G, cellID);
    scalar dfdq = (fxplusq - fxminusq)/(2*dx);
    // ToDo: check analytical derivative

    // New Guess for x
    residual = (fdeq/dfdq);
    deq -= residual;
    residual /= maxMagDEpsilon;

    // fdeq will go to zero at convergence
    fdeq = gfun(ebart, deq, qe, G, cellID);

    if (i > MaxNewtonIter_-2)
      {
        Warning
            << "Aravas plasticity not converging, fx is " << mag(fdeq) <<endl;
      }
      }
    // while ( (mag(fdeq) > LoopTol_) && ++i < MaxNewtonIter_);
    while
        (
            (mag(residual) > LoopTol_)
            && ++i < MaxNewtonIter_
            ); // use change in deq instead

    // update yield stress
    scalar q = qfun(qe, deq, G);
    s0 = s0fun(ebart, cellID);
    scalar ebar = ebarfun(ebart, q, deq, s0);
    s0 = s0fun(ebar, cellID);
  }

void aravasMises::updateYieldStress()
{
  Info << nl << "Updating the yield stress" << endl;
  sigmaY_ += DSigmaY_;
  sigmaY_.correctBoundaryConditions();

  // count cells actively yielding
  int numCellsYielding = 0;
  forAll(activeYield_.internalField(), celli)
    {
      if (DEpsilonPEq_.internalField()[celli] > SMALL)
    {
      activeYield_.internalField()[celli] = 1.0;
      numCellsYielding++;
    }
      else
    {
      activeYield_.internalField()[celli] = 0.0;
    }
    }
  reduce(numCellsYielding, sumOp<int>());
  forAll(activeYield_.boundaryField(), patchi)
    {
      if (!activeYield_.boundaryField()[patchi].coupled())
    {
      forAll(activeYield_.boundaryField()[patchi], facei)
        {
          if (DEpsilonPEq_.boundaryField()[patchi][facei] > SMALL)
        {
          activeYield_.boundaryField()[patchi][facei] = 1.0;
        }
          else
        {
          activeYield_.boundaryField()[patchi][facei] = 0.0;
        }
        }
    }
    }
  Info << "\t" << numCellsYielding << " cells are actively yielding"
       << nl << endl;


  Info << "\tMax DEpsilonPEq is " << gMax(DEpsilonPEq_) << nl << endl;
}

} // end of namespace
// ************************************************************************* //
