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

#include "yamadaMises.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "constitutiveModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(yamadaMises, 0);
    addToRunTimeSelectionTable(plasticityStressReturn, yamadaMises, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
yamadaMises::yamadaMises
(
 const word& name,
 constitutiveModel& constitutiveModel
)
:
  plasticityStressReturn(name, constitutiveModel),
  constitutiveModel_(constitutiveModel),
  sigmaY_(constitutiveModel.sigmaY()),
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
  epsilonPEq_
  (
   IOobject
   (
    "epsilonPEq",
    sigmaY_.time().timeName(),
    sigmaY_.db(),
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
    ),
   sigmaY_.mesh(),
   dimensionedScalar("0", dimless, 0.0)
   ),
  beta_
  (
   IOobject
   (
    "beta",
    sigmaY_.time().timeName(),
    sigmaY_.db(),
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
    ),
   sigmaY_.mesh(),
   dimensionedScalar("0", dimless, 0)
   )
{
  Info << "Creating yamadaMises plasticity stress return method" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

yamadaMises::~yamadaMises()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void yamadaMises::correct()
{
  //Info << "\tCorrecting plasticity model ... " << flush;

  // we will lookup the strain increment for the solver
  // instead of calculating it as it has already been calculated

  // for under-relaxation
  DEpsilonP_.storePrevIter();

  // mesh
  const fvMesh& mesh = sigmaY_.mesh();

  // lookup mechanical properties
  const volScalarField& mu =
    mesh.objectRegistry::lookupObject<volScalarField>("mu");
  const volScalarField& lambda =
    mesh.objectRegistry::lookupObject<volScalarField>("lambda");

  // Lookup current strain increment from the solver
  const volSymmTensorField& DEpsilon =
    mesh.objectRegistry::lookupObject<volSymmTensorField>("DEpsilon");

  // old total strain
  const volSymmTensorField& epsilon =
    mesh.objectRegistry::lookupObject<volSymmTensorField>("epsilon");

  // Equivalent strain
    const volScalarField epsilonEq =
        sqrt((2.0/3.0)*magSqr(dev(epsilon + DEpsilon)))
      + dimensionedScalar("SMALL", dimless, SMALL);


    // Current yield stress
    scalarField& sigmaYI = sigmaY_.internalField();
    // Initial yield stress
    scalarField initialSigmaYI = constitutiveModel_.sigmaY()().internalField();

    volScalarField epsilonEqCorr = epsilonEq;

    // Update mu and lambda - for nonlinear elastic
    //mu= mu(epsilonEqCorr);
    //lambda = lambda(epsilonEqCorr);

    const volScalarField DEpsilonEq =
        sqrt((2.0/3.0)*magSqr(dev(epsilon + DEpsilon)))
      - sqrt((2.0/3.0)*magSqr(dev(epsilon)))
      + dimensionedScalar("SMALL", dimless, SMALL);

    const volSymmTensorField DSigma =
        2*mu*(DEpsilon - DEpsilonP_) + I*(lambda*tr(DEpsilon));
    //const volSymmTensorField& oldSigma = sigma();
    const volSymmTensorField& oldSigma =
      mesh.objectRegistry::lookupObject<volSymmTensorField>("sigma");
    volScalarField oldSigmaEq("oldSigmaEq", sqrt(1.5*magSqr(dev(oldSigma))));

    const volSymmTensorField newSigma = oldSigma + DSigma;
    const volScalarField sigmaEq
      (
       "sigmaEq",
       sqrt(1.5*magSqr(dev(newSigma)))
       + dimensionedScalar("SMALL", dimPressure, SMALL)
       );
    const volSymmTensorField devSigma = dev(newSigma);
    const volSymmTensorField DSigmaE = DSigma + 2*mu*DEpsilonP_;
    const volScalarField sigmaEqE  = sqrt(1.5*magSqr(dev(oldSigma + DSigmaE)));
    const volScalarField DSigmaEqE = sqrt(1.5*magSqr(dev(DSigmaE)));

    // old total plastic strain
    const volSymmTensorField& oldEpsilonP =
      mesh.objectRegistry::lookupObject<volSymmTensorField>("epsilonP");
    // plastic equivalent strain
    const volScalarField epsilonPEq =
        sqrt((2.0/3.0)*magSqr(dev(oldEpsilonP+DEpsilonP_)));
    const volScalarField oldEpsilonPEq
        ("oldEpsilonPEq", sqrt((2.0/3.0)*magSqr(dev(oldEpsilonP))));

    // Get current plastic modulus which is, in general,
    // a function of epsilonPEq (plastic equivalent strain)
    // Ep is constant for perfect/linear plastic
    // nonLinearElasticPlastic calculates Ep using
    // sigmaEq
    //volScalarField Ep_ = Ep(sigmaEq);
    volScalarField Ep_ = constitutiveModel_.Ep(epsilonPEq);
    //volScalarField oldEp = Ep(oldEpsilonPEq);
    //volScalarField newEp = Ep(epsilonPEq);
    //scalar weightEp = 0.75;
    //volScalarField Ep_ = weightEp*newEp + (1.0-weightEp)*oldEp;

    // Update internal beta
    const scalarField& muI = mu.internalField();
    const scalarField& lambdaI = lambda.internalField();
    const symmTensorField& DEpsilonI = DEpsilon.internalField();
    const scalarField& DEpsilonEqI  = DEpsilonEq.internalField();
    const symmTensorField& oldSigmaI = oldSigma.internalField();
    const scalarField& oldSigmaEqI = oldSigmaEq.internalField();
    const symmTensorField& devSigmaI = devSigma.internalField();
    const symmTensorField& DSigmaEI = DSigmaE.internalField();
    const scalarField& sigmaEqEI  = sigmaEqE.internalField();
    //const scalarField& sigmaEqI  = sigmaEq.internalField();
    const scalarField& DSigmaEqEI = DSigmaEqE.internalField();
    const scalarField& oldBetaI = beta_.oldTime().internalField();

    scalarField& betaI = beta_.internalField();

    forAll (betaI, cellI)
    {
        tensor curDEpsEPred = tensor::zero;

        if ( (DEpsilonEqI[cellI] >= 0) && (oldBetaI[cellI] > SMALL) )
        {
            // philipc: this means that once beta=1 then it can never change
            // back to zero as DEpsilonEq will probably always be finite
            // I will add a check here if the cell has unloaded from the
            // yield surface
            //betaI[cellI] = 1.0;
            curDEpsEPred = tensor::zero;

        // philipc
        if (sigmaEqEI[cellI] < (sigmaYI[cellI]-SMALL))
          {
        betaI[cellI] = 0.0;
          }
        else
          {
        betaI[cellI] = 1.0;
          }
        }
        else
        {
            betaI[cellI] = 0.0;
            curDEpsEPred = DEpsilonI[cellI];

            if
            (
                (DEpsilonEqI[cellI] >= 0)
        && (sigmaEqEI[cellI] > (sigmaYI[cellI]-SMALL))
            )
            {
                scalar C = sqr(oldSigmaEqI[cellI]) - sqr(sigmaYI[cellI]);
                scalar B = 3.0*(dev(oldSigmaI[cellI]) && dev(DSigmaEI[cellI]));
                scalar A = sqr(DSigmaEqEI[cellI]);

        scalar alpha = (-B + ::sqrt(mag(B*B - 4*A*C)))/(2*A + SMALL);
                //scalar alpha = (-B + ::sqrt((B*B - 4*A*C)))/(2*A + SMALL);
                curDEpsEPred =
                    alpha/(2.0*muI[cellI] + SMALL)
                   *(
                        DSigmaEI[cellI]
                      - (lambdaI[cellI]
                         /(2*muI[cellI] + 3*lambdaI[cellI] + SMALL))
                       *tr(DSigmaEI[cellI])*I
                    );

                betaI[cellI] =
                    1.0
                  - (devSigmaI[cellI] && curDEpsEPred)
           /((devSigmaI[cellI] && DEpsilonI[cellI]) + SMALL);
            }
        }

        betaI[cellI] = max(betaI[cellI], 0.0);
        betaI[cellI] = min(betaI[cellI], 1.0);
    }

    // Update beta at boundary
    forAll(beta_.boundaryField(), patchI)
      {
        if (!beta_.boundaryField()[patchI].coupled())
      {
        const scalarField& muPatch = mu.boundaryField()[patchI];
        const scalarField& lambdaPatch = lambda.boundaryField()[patchI];
        const scalarField& sigmaYPatch = sigmaY_.boundaryField()[patchI];
        const symmTensorField& DEpsilonPatch =
          DEpsilon.boundaryField()[patchI];
        const scalarField DEpsilonEqPatch =
          DEpsilonEq.boundaryField()[patchI];
        const symmTensorField& oldSigmaPatch =
          oldSigma.boundaryField()[patchI];
        const scalarField& oldSigmaEqPatch =
          oldSigmaEq.boundaryField()[patchI];
        const symmTensorField& devSigmaPatch =
          devSigma.boundaryField()[patchI];
        const symmTensorField& DSigmaEPatch = DSigmaE.boundaryField()[patchI];
        const scalarField& sigmaEqEPatch = sigmaEqE.boundaryField()[patchI];
        //const scalarField& sigmaEqPatch = sigmaEq.boundaryField()[patchI];
        const scalarField& DSigmaEqEPatch = DSigmaEqE.boundaryField()[patchI];
        const scalarField& oldBetaPatch =
          beta_.oldTime().boundaryField()[patchI];
        scalarField& betaPatch = beta_.boundaryField()[patchI];

        forAll(betaPatch, faceI)
          {
        tensor curDEpsEPred = tensor::zero;
        if
          (
           (DEpsilonEqPatch[faceI] >= 0)
           && (oldBetaPatch[faceI] > SMALL)
            )
          {
            //betaPatch[faceI] = 1;
            curDEpsEPred = tensor::zero;

            // philipc
            if (sigmaEqEPatch[faceI] < (sigmaYPatch[faceI]-SMALL))
              {
            betaPatch[faceI] = 0.0;
              }
            else
              {
            betaPatch[faceI] = 1;
              }
          }
        else
          {
            betaPatch[faceI] = 0;
            curDEpsEPred = DEpsilonPatch[faceI];

            if
              (
               (DEpsilonEqPatch[faceI] >= 0)
               && (sigmaEqEPatch[faceI] > (sigmaYPatch[faceI]-SMALL))
            )
              {
            scalar C =
              sqr(oldSigmaEqPatch[faceI])
              - sqr(sigmaYPatch[faceI]);
            scalar B =
              3.0
              *(
                            dev(oldSigmaPatch[faceI])
                && dev(DSigmaEPatch[faceI])
                );
            scalar A = sqr(DSigmaEqEPatch[faceI]);
            scalar alpha = (-B + ::sqrt(mag(B*B-4*A*C)))/(2*A + SMALL);
            //scalar alpha = (-B + ::sqrt((B*B-4*A*C)))/(2*A + SMALL);

            curDEpsEPred =
              alpha/(2.0*muPatch[faceI] + SMALL)
              *(
                            DSigmaEPatch[faceI]
                - (
                   lambdaPatch[faceI]
                               /
                   (
                       2*muPatch[faceI]
                       + 3*lambdaPatch[faceI] + SMALL
                       )
                   )
                *tr(DSigmaEPatch[faceI])*I
                );

            betaPatch[faceI] =
              1.0
              - (devSigmaPatch[faceI] && curDEpsEPred)
              /((devSigmaPatch[faceI] && DEpsilonPatch[faceI]) + SMALL);
              }
          }

        betaPatch[faceI] = max(betaPatch[faceI], 0.0);
        betaPatch[faceI] = min(betaPatch[faceI], 1.0);
          }
      }
      }

    beta_.correctBoundaryConditions();

    // Update plastic strain increment
    DEpsilonP_ =
      4.5*beta_*mu*(devSigma && DEpsilon)*devSigma
      /(
    (Ep_ + 3*mu)*sqr(sigmaEq)
    + dimensioned<scalar>
    (
     "SMALL",
     mu.dimensions()*sigmaEq.dimensions()*sigmaEq.dimensions(),
     SMALL
     )
        );

    // DEpsilonP_ = rf*newDEpsilonP + (1.0 - rf)*DEpsilonP_;

    DEpsilonP_.relax();
    DEpsilonP_.correctBoundaryConditions();

    // forAll(DEpsilonP_, celli)
    //   {
    //      if (mag(DEpsilonP_[celli]) > SMALL && beta_[celli] < SMALL)
    //        {
    //          DEpsilonP_[celli] = symmTensor::zero;
    //        }
    //   }
    // forAll(DEpsilonP_.boundaryField(), patchi)
    //   {
    //      forAll(DEpsilonP_.boundaryField()[patchi], facei)
    //        {
    //          if (mag(DEpsilonP_.boundaryField()[patchi][facei]) > SMALL
    //             && beta_.boundaryField()[patchi][facei] < SMALL)
    //            {
    //          DEpsilonP_.boundaryField()[patchi][facei] = symmTensor::zero;
    //            }
    //        }
    //   }

    //Info << "done" << endl;
}

  void yamadaMises::updateYieldStress()
  {
    Info << nl << "Updating the yield stress" << endl;

    // update total equivalent plastic strain
    epsilonPEq_ += DEpsilonPEq_;

    // mesh
    const fvMesh& mesh = sigmaY_.mesh();

    //const volSymmTensorField& newSigma = sigma();
    //const volSymmTensorField& newSigma =
    //mesh.objectRegistry::lookupObject<volSymmTensorField>("sigma");
    //const volScalarField sigmaEq = sqrt(1.5*magSqr(dev(newSigma)));

    // old total plastic strain
    //const volSymmTensorField& oldEpsilonP =
    //mesh.objectRegistry::lookupObject<volSymmTensorField>("epsilonP");
    // plastic equivalent strain
    //const volScalarField epsilonPEq =
    //sqrt((2.0/3.0)*magSqr(dev(oldEpsilonP+DEpsilonP_)));
    //const volScalarField oldEpsilonPEq
    //("oldEpsilonPEq",sqrt((2.0/3.0)*magSqr(dev(oldEpsilonP))));

    // Get current plastic modulus which is, in general,
    // a function of epsilonPEq (plastic equivalent strain)
    // Ep is constant for perfect/linear plastic
    // nonLinearElasticPlastic calculates Ep using
    // sigmaEq
    //volScalarField Ep_ = Ep(sigmaEq);
    //const volSymmTensorField& DSigma =
    //mesh.objectRegistry::lookupObject<volSymmTensorField>("DSigma");
    // oldSigmaEq is looked up by rheologyLaw
    //volScalarField oldSigmaEq
    //("oldSigmaEq", sqrt(1.5*magSqr(dev(newSigma-DSigma))));
    //volScalarField Ep_ = constitutiveModel_.Ep(epsilonPEq);
    // volScalarField oldEp_ = Ep(oldEpsilonPEq);
    // volScalarField newEp_ = Ep(epsilonPEq);
    // volScalarField Ep_ = 0.5*oldEp_ + 0.5*newEp_;

    //const scalarField& EpI = Ep_.internalField();
    //const scalarField& sigmaEqI = sigmaEq.internalField();

    const scalarField& betaI = beta_.internalField();
    const scalarField& epsilonPEqI = epsilonPEq_.internalField();
    scalarField& sigmaYI = sigmaY_.internalField();

    label numCellsUpdated = 0;

    forAll(sigmaYI, cellI)
    {
      // philipc: sigmaY should be updated based on current
      // total plastic equivalent strain
      if (betaI[cellI] > SMALL)
    {
      sigmaYI[cellI] = constitutiveModel_.sigmaY(epsilonPEqI[cellI], cellI);
    }

        // if (EpI[cellI] != 0)
        // {
        //     if ( sigmaEqI[cellI] > sigmaYI[cellI] )
        //     {
        //         sigmaYI[cellI] = sigmaEqI[cellI];
    //      numCellsUpdated++;
    // }
    //  }
    }

    forAll(sigmaY_.boundaryField(), patchI)
      {
        if (!sigmaY_.boundaryField()[patchI].coupled())
      {
        //const scalarField& EpPatch = Ep_.boundaryField()[patchI];
        //const scalarField& sigmaEqPatch = sigmaEq.boundaryField()[patchI];
        const scalarField& betaPatch = beta_.boundaryField()[patchI];
        const scalarField& epsilonPEqPatch =
            epsilonPEq_.boundaryField()[patchI];
        scalarField& sigmaYPatch = sigmaY_.boundaryField()[patchI];
        const labelList& faceCells = mesh.boundaryMesh()[patchI].faceCells();

        forAll(sigmaYPatch, faceI)
          {
        if (betaPatch[faceI] > SMALL)
          {
            // ID is for multiMaterial
            sigmaYPatch[faceI] =
                constitutiveModel_.sigmaY
                (epsilonPEqPatch[faceI], faceCells[faceI]);
          }

        // if (EpPatch[faceI] != 0)
        // {
        //     if (sigmaEqPatch[faceI] > sigmaYPatch[faceI])
        //     {
        //         sigmaYPatch[faceI] = sigmaEqPatch[faceI];

        //         //Info << "Boundary cell " << patchI << " " << faceI
        //          //  << " Yield stress updated to Sy= "
        //          //  << sigmaEqPatch[faceI] * 1.0E-06 << " MPa"
        //          //  << endl;
        //     }
        // }
          }
      }
      }

    Info << "\t" << numCellsUpdated << " cells have been updated"
     << endl;


    // count cells actively yielding
    int numCellsYielding = 0;
    forAll(beta_.internalField(), celli)
      {
    if (beta_.internalField()[celli] > SMALL)
      numCellsYielding++;
      }
    reduce(numCellsYielding, sumOp<int>());
    Info << "\t" << numCellsYielding << " cells are actively yielding"
     << endl;


    // we will also set DEpsilonPEq which is just for visualisation
    DEpsilonPEq_ = sqrt((2.0/3.0)*magSqr(dev(DEpsilonP_)));

    Info << "\tMax DEpsilonPEq is " << gMax(DEpsilonPEq_) << nl << endl;

    //Info << "done" << endl;
}

} // end of namespace
// ************************************************************************* //
