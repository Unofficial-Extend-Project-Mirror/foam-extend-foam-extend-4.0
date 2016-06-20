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

#include "newAravasMises.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "calculatedFvPatchFields.H"
#include "constitutiveModel.H"
#include "solidSolver.H"
// #include "pRveUnsIncrTotalLagrangianSolid.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(newAravasMises, 0);
    addToRunTimeSelectionTable
    (
        plasticityStressReturn,
        newAravasMises,
        dictionary
    );


    // Static variables

        // Tolerance for Newton loop
        scalar newAravasMises::LoopTol_ = 1e-8;

        // Maximum number of iterations for Newton loop
        label newAravasMises::MaxNewtonIter_ = 100;

        // finiteDiff is the delta for finite difference differentiation
        scalar newAravasMises::finiteDiff_ = 0.25e-6;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


scalar newAravasMises::qfun
(
    const scalar qe,
    const scalar deq,
    const scalar G
    ) const
{
    return (qe - 3*G*deq);
}

// yield stress
scalar newAravasMises::s0fun (const scalar ebar, const label cellID) const
{
    scalar ebarNew = ebar;
    if (ebarNew < SMALL)
    {
        ebarNew = SMALL;
    }

    return
        constitutiveModel_.sigmaY
        (
            ebarNew,
            cellID
        ); //strainYieldStressXY(ebarNew);
}

scalar newAravasMises::ebarfun
(
    const scalar ebart,
    const scalar q,
    const scalar deq,
    const scalar s0
) const
{
    // Optimise integration for large steps
    return (ebart + (q*deq/s0));

    // ZT ?
//     return (ebart + deq);
}

scalar newAravasMises::gfun
(
    const scalar ebart,
    const scalar deq,
    const scalar qe,
    const scalar G,
    const label cellID
) const
{
    const scalar q = qfun(qe, deq, G);

    scalar s0 = s0fun(ebart + deq, cellID);

//     scalar s0 = s0fun(ebart, cellID);
//     const scalar ebar = ebarfun(ebart, q, deq, s0);
//     s0 = s0fun(ebar, cellID);

    return (Foam::pow(q/s0, 2) - 1.0);
}

void newAravasMises::aravasNewtonLoop
(
    scalar& deq,
    scalar& s0,
    //const scalar s00,
    const scalar G,
    const scalar ebart,
    const scalar qe,
    const label cellID,
    const scalar maxMagDEpsilon
) const
{
    // Loop to determine DEpsilonPEq
    // using Newtion's method

    /*
      J2 Mises Yield function:
      f = (q/s0)^2 - 1.0 = 0.0

      Substituting "q = qe - 3*G*deq" into f:
      f = ((qe - 3*G*deq)/s0)^2 - 1.0 = 0.0
      where
      s0 is the current yield stress and is a function of deq;
      deq is the current increment of plastic equivalent strain;
      q is the equivalent stress where "q = qe - 3*G*deq";
      qe is the trial elastic equivalent stress;
      G is the shear modulus.

      It is not straight-forward to calculate df/d(deq) analytically,
      so we will numerically determine the df/d(deq) using finite differences.

      Notes on numerical differentiation, from Hauser.
      Although the Second Order derivative seems better, the First Order
      approximation is actually more efficient when used in Newton's method
      as less function evaluations are required and the converged value is not
      affected by the accuracy of the derivative approximation.
      Hauser also gives considerable emphasis on the choice of dx and suggest
      dx = epsilon*x where epsilon = 1e-6;

      First Order
      Require 2 function evaluations, but one of those corresponds
      with function evaluation in Newton's method
      f'(x) = (f(x) - f(x - dx))/(dx)

      Second Order
      Require 2 function evaluations
      f'(x) = (f(x + dx) - f(x - dx))/(2*dx)

      Fourth Order
      Require 4 function evaluations
      f'(x) = (-f(x + 2*dx) + 8*f(x + dx) - 8*f(x - dx) + f(x + 2*dx))/(12*dx)
    */

    int i = 0;
    scalar fdeq = gfun(ebart, deq, qe, G, cellID);
    scalar residual = 1.0;
    do
    {
        // Numerically calculate parial derivatives

        // First Order
        // Hauser 2009 suggested First Order is more efficient for Newton's
        // Method as only two function evaluaitons are required
        scalar fxplusq  = gfun(ebart, deq + finiteDiff_, qe, G, cellID);
        scalar dfdq = (fxplusq - fdeq)/finiteDiff_;

        // Second Order
        //scalar fxplusq  = gfun(ebart, deq+finiteDiff_, qe, G, cellID);
        //scalar fxminusq  = gfun(ebart, deq-finiteDiff_, qe, G, cellID);
        //scalar dfdq = (fxplusq - fxminusq)/(2*finiteDiff_);

        // New Guess for x
        residual = (fdeq/dfdq);
        deq -= residual;
        residual /= maxMagDEpsilon; // Normalise wrt strain increment

        // fdeq will go to zero at convergence
        fdeq = gfun(ebart, deq, qe, G, cellID);

        if (i == MaxNewtonIter_)
        {
            Warning << "Aravas plasticity not converging, fx is "
                << mag(fdeq) <<endl;
        }
    }
    // while ( (mag(fdeq) > LoopTol_) && ++i < MaxNewtonIter_);
    while( (mag(residual) > LoopTol_) && ++i < MaxNewtonIter_);

    // update yield stress

//     scalar q = qfun(qe, deq, G);
//     s0 = s0fun(ebart, cellID);
//     scalar ebar = ebarfun(ebart, q, deq, s0);
//     s0 = s0fun(ebar, cellID);

    s0 = s0fun(ebart+deq, cellID);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
newAravasMises::newAravasMises
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
    sigmaYf_
    (
        IOobject
        (
            "sigmaYf",
            constitutiveModel.sigma().time().timeName(),
            constitutiveModel.sigma().db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(constitutiveModel.sigmaY())
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
    DSigmaYf_
    (
        IOobject
        (
            "DSigmaYf",
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
    epsilonP_
    (
        IOobject
        (
            "epsilonP",
            sigmaY_.time().timeName(),
            sigmaY_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sigmaY_.mesh(),
        dimensionedSymmTensor("0", dimless, symmTensor::zero)
    ),
    DEpsilonPf_
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
    DEpsilonPEqf_
    (
        IOobject
        (
            "DEpsilonPEqf",
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
    epsilonPEqf_
    (
        IOobject
        (
            "epsilonPEqf",
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
    ),
    plasticNf_
    (
        IOobject
        (
            "plasticNf",
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

newAravasMises::~newAravasMises()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void newAravasMises::correct()
{
    // Fields for integration of plasticity model

    // mesh
    const fvMesh& mesh = sigmaY_.mesh();

    symmTensor avgDEpsilon = symmTensor::zero;
    symmTensor avgEpsilon = symmTensor::zero;

    // // Looking up solid solver
    // const solidSolver& stress =
    //     mesh.lookupObject<solidSolver>
    //     (
    //         "solidProperties"
    //     );

//     if (isA<solidSolvers::pRveUnsIncrTotalLagrangianSolid>(stress))
//     {
//         // Get average deformation gradient from the solver
//         const solidSolvers::pRveUnsIncrTotalLagrangianSolid& pRveStress =
//             refCast<const solidSolvers::pRveUnsIncrTotalLagrangianSolid>
//             (
//                 stress
//             );

//         avgDEpsilon = pRveStress.avgDEpsilon();
//         avgEpsilon = pRveStress.avgEpsilon();
//     }

    // Cell plastic strain
    if (true)
    {
        // elastic properties
        const volScalarField& mu =
            mesh.objectRegistry::lookupObject<volScalarField>("mu");
        const volScalarField& lambda =
            mesh.objectRegistry::lookupObject<volScalarField>("lambda");

        // Old plastic strain
        const volSymmTensorField& oldEpsilonP =
            mesh.objectRegistry::lookupObject<volSymmTensorField>("epsilonP");

        // Old plastic equivalent strain
        const volScalarField oldEpsilonPEq =
            sqrt((2.0/3.0)*magSqr(dev(oldEpsilonP)));

        // Current strain increment - should be recalculated within
        // momentum loop before rheology.correct()
        const volSymmTensorField& DEpsilon =
            mesh.objectRegistry::lookupObject<volSymmTensorField>("DEpsilon");

        // Old total strain
        const volSymmTensorField& oldEpsilon =
            mesh.objectRegistry::lookupObject<volSymmTensorField>("epsilon");

        // New total strain
        const volSymmTensorField newEpsilon =
            oldEpsilon + avgEpsilon + DEpsilon + avgDEpsilon;

        // Magnitude of current strain increment
        volScalarField magDEpsilon = Foam::magSqr(DEpsilon + avgDEpsilon);

        // Calculate Old Elastic Strain
        volSymmTensorField oldEpsilonElastic =
            oldEpsilon + avgEpsilon - oldEpsilonP;

        // Calculate old Elastic Stress
        volSymmTensorField oldSigmaElastic =
            2.0*mu*oldEpsilonElastic + lambda*tr(oldEpsilonElastic)*I;

        // Calculate old Elastic Stress deviatoric component
        volSymmTensorField oldSigmaElasticDev = dev(oldSigmaElastic);

        // Calculate old equivalent elastic stress
        volScalarField oldSigmaEqElastic =
            sqrt((3.0/2.0)*Foam::magSqr(oldSigmaElasticDev));

        // Calculate new Elastic Strain - this is a trial elastic strain
        volSymmTensorField newEpsilonElastic = newEpsilon - oldEpsilonP;

        // Calculate new Elastic Stress
        volSymmTensorField newSigmaElastic =
            2.0*mu*newEpsilonElastic + lambda*tr(newEpsilonElastic)*I;

        // Calculate change in elastic Stress
        volSymmTensorField DSigmaElastic = newSigmaElastic - oldSigmaElastic;

        // Calculate new elastic Deviatoric
        volSymmTensorField newSigmaDevElastic = dev(newSigmaElastic);

        // Calculate new elastic equivalent stress
        volScalarField newSigmaEqElastic =
            sqrt((3.0/2.0)*magSqr(newSigmaDevElastic));

        // Store previous iteration for under-relaxation
        DEpsilonP_.storePrevIter();

        // Normalise residual in Newton method with respect to DEpsilon
        const scalar maxMagDEpsilon =
            max(gMax(mag(DEpsilon.internalField())), SMALL);

        // Calculate beta field (fully elastic-plastic fraction)
        //# include "newAravasMisesUpdateBeta.H"

        // Calculate DEpsilonPEq_ and plasticN_
        forAll (DEpsilonP_, cellI)
        {
            // Calculate return direction plasticN

            if(newSigmaEqElastic[cellI] > SMALL)
            {
                plasticN_[cellI] =
                    3.0*newSigmaDevElastic[cellI]
                   /(2.0*newSigmaEqElastic[cellI]);
            }

            // J2 yield function
            scalar fy =
                Foam::pow
                (
                    newSigmaEqElastic[cellI]/sigmaY_[cellI],
                    2.0
                )
              - 1.0;

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
                (
                    DEpsilonPEq_[cellI],
                    s0,
                    //s00,
                    G,
                    ebart,
                    qe,
                    cellI,
                    maxMagDEpsilon
                );

                // Update equivalent plastic strain and
                // increment of yield stress
                DSigmaY_[cellI] = s0 - s00;
            }
        }

        forAll(plasticN_.boundaryField(), patchI)
        {
            if (!plasticN_.boundaryField()[patchI].coupled())
            {
                const labelList& faceCells =
                    mesh.boundary()[patchI].faceCells();

                forAll(plasticN_.boundaryField()[patchI], faceI)
                {
                    // Calculate direction plasticN
                    if
                    (
                        newSigmaEqElastic.boundaryField()[patchI][faceI]
                      > SMALL
                    )
                    {
                        plasticN_.boundaryField()[patchI][faceI] =
                        (
                            3.0
                           *newSigmaDevElastic.boundaryField()[patchI][faceI]
                           /(
                                2.0
                               *newSigmaEqElastic.boundaryField()
                                [patchI][faceI]
                            )
                        );
                    }

                    // J2 yield function
                    scalar fy =
                        Foam::pow
                        (
                            newSigmaEqElastic.boundaryField()[patchI][faceI]
                           /sigmaY_.boundaryField()[patchI][faceI],
                            2.0
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

                        // total equivalent strain matrix
                        // (t is start of time-step)
                        const scalar ebart =
                            oldEpsilonPEq.boundaryField()[patchI][faceI];
                        // yield stress at start of time-step
                        const scalar s00 =
                            sigmaY_.boundaryField()[patchI][faceI];
                        scalar s0 = 0.0; // updated in loop below
                        const scalar qe =
                            newSigmaEqElastic.boundaryField()[patchI][faceI];
                        const scalar G = mu.boundaryField()[patchI][faceI];

                        // Calculate deq and s0
                        aravasNewtonLoop
                        (
                            DEpsilonPEq_.boundaryField()[patchI][faceI],
                            s0,
                            //s00,
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

        DSigmaY_.correctBoundaryConditions();
        DEpsilonPEq_.correctBoundaryConditions();
        plasticN_.correctBoundaryConditions();

        // Update plastic strain increment
        DEpsilonP_ = DEpsilonPEq_*plasticN_;

//         forAll(DEpsilonP_.boundaryField(), patchI)
//         {
//             if
//             (
//                 DEpsilonP_.boundaryField()[patchI].type()
//              == calculatedFvPatchSymmTensorField::typeName
//             )
//             {
//                 DEpsilonP_.boundaryField()[patchI] =
//                     DEpsilonP_.boundaryField()[patchI].patchInternalField();
//             }
//         }

        // Relax plastic strain increment to help convergence
        DEpsilonP_.relax();
    }

    // Face plastic strain
    if (false)
    {
        // elastic properties
        const surfaceScalarField& muf =
            mesh.objectRegistry::lookupObject<surfaceScalarField>("muf");
        const surfaceScalarField& lambdaf =
            mesh.objectRegistry::lookupObject<surfaceScalarField>("lambdaf");

        // Old plastic strain
        const surfaceSymmTensorField& oldEpsilonPf =
            mesh.objectRegistry::lookupObject<surfaceSymmTensorField>
            (
                "epsilonPf"
            );

        // Old plastic equivalent strain
        const surfaceScalarField oldEpsilonPEqf =
            sqrt((2.0/3.0)*magSqr(dev(oldEpsilonPf)));

        // Current strain increment - should be recalculated within
        // momentum loop before rheology.correct()
        const surfaceSymmTensorField& DEpsilonf =
            mesh.objectRegistry::lookupObject<surfaceSymmTensorField>
            (
                "DEpsilonf"
            );

        // Old total strain
        const surfaceSymmTensorField& oldEpsilonf =
            mesh.objectRegistry::lookupObject<surfaceSymmTensorField>
            (
                "epsilonf"
            );

        // New total strain
        const surfaceSymmTensorField newEpsilonf =
            oldEpsilonf + avgEpsilon + DEpsilonf + avgDEpsilon;

        // Magnitude of current strain increment
        surfaceScalarField magDEpsilonf =
            Foam::magSqr(DEpsilonf + avgDEpsilon);

        // Calculate Old Elastic Strain
        surfaceSymmTensorField oldEpsilonElasticf =
            oldEpsilonf + avgEpsilon - oldEpsilonPf;

        // Calculate old Elastic Stress
        surfaceSymmTensorField oldSigmaElasticf =
            2.0*muf*oldEpsilonElasticf + lambdaf*tr(oldEpsilonElasticf)*I;

        // Calculate old Elastic Stress deviatoric component
        surfaceSymmTensorField oldSigmaDevf = dev(oldSigmaElasticf);

        // Calculate old equivalent stress
        surfaceScalarField oldSigmaEqElasticf =
            sqrt((3.0/2.0)*Foam::magSqr(oldSigmaDevf));

        // Calculate new Elastic Strain - this is a trial elastic stress
        surfaceSymmTensorField newEpsilonElasticf =
            newEpsilonf  - oldEpsilonPf;

        // Calculate new Elastic Stress
        surfaceSymmTensorField newSigmaElasticf =
            2.0*muf*newEpsilonElasticf + lambdaf*tr(newEpsilonElasticf)*I;

        // Calculate change in elastic Stress
        surfaceSymmTensorField DSigmaElasticf =
            newSigmaElasticf - oldSigmaElasticf;

        // Calculate new elastic Deviatoric
        surfaceSymmTensorField newSigmaDevElasticf = dev(newSigmaElasticf);

        // Calculate new elastic equivalent stress
        surfaceScalarField newSigmaEqElasticf =
            sqrt((3.0/2.0)*magSqr(newSigmaDevElasticf));

        // Store previous iteration for under-relaxation
        DEpsilonPf_.storePrevIter();

        // Normalise residual in Newton method with respect to DEpsilon
        const scalar maxMagDEpsilon =
            max(gMax(mag(DEpsilonf.internalField())), SMALL);

        // Calculate DEpsilonPEq_ and plasticN_
        forAll (DEpsilonPf_, faceI)
        {
            // Calculate return direction plasticN

            if (newSigmaEqElasticf[faceI] > SMALL)
            {
                plasticNf_[faceI] =
                    3.0*newSigmaDevElasticf[faceI]
                   /(2.0*newSigmaEqElasticf[faceI]);
            }

            // J2 yield function
            scalar fy =
                Foam::pow
                (
                    newSigmaEqElasticf[faceI]/sigmaYf_[faceI],
                    2.0
                )
              - 1.0;

            // Calculate DEpsilonPEq
            if (fy < SMALL)
            {
                // elastic
                DSigmaYf_[faceI] = 0.0;
                DEpsilonPEqf_[faceI] = 0.0;
            }
            else
            {
                // yielding

                // total equivalent strain matrix (t is start of time-step)
                const scalar ebart = oldEpsilonPEqf[faceI];

                // yield stress at start of time-step
                const scalar s00 = sigmaYf_[faceI];
                scalar s0 = 0.0; // updated in loop below
                const scalar qe = newSigmaEqElasticf[faceI];
                const scalar G = muf[faceI];

                // This Aravas loop calculates deq using Newtons's method
                // Calculate deq and s0
                aravasNewtonLoop
                (
                    DEpsilonPEqf_[faceI],
                    s0,
                    //s00,
                    G,
                    ebart,
                    qe,
                    faceI,
                    maxMagDEpsilon
                );

                // Update equivalent plastic strain and
                // increment of yield stress
                DSigmaYf_[faceI] = s0 - s00;
            }
        }

        forAll(plasticNf_.boundaryField(), patchI)
        {
//             if (!plasticNf_.boundaryField()[patchI].coupled())
            {
                const labelList& faceCells =
                    mesh.boundary()[patchI].faceCells();

                forAll(plasticNf_.boundaryField()[patchI], faceI)
                {
                    // Calculate direction plasticN
                    if
                    (
                        newSigmaEqElasticf.boundaryField()[patchI][faceI]
                      > SMALL
                    )
                    {
                        plasticNf_.boundaryField()[patchI][faceI] =
                        (
                            3.0
                           *newSigmaDevElasticf.boundaryField()[patchI][faceI]
                           /(
                                2.0
                               *newSigmaEqElasticf.boundaryField()
                                [patchI][faceI]
                            )
                        );
                    }

                    // J2 yield function
                    scalar fy =
                        Foam::pow
                        (
                            newSigmaEqElasticf.boundaryField()[patchI][faceI]
                           /sigmaYf_.boundaryField()[patchI][faceI],
                            2.0
                        )
                      - 1.0;

                    // Calculate DEpsilonPEq
                    if (fy < SMALL)
                    {
                        // elasticity
                        DSigmaYf_.boundaryField()[patchI][faceI] = 0.0;
                        DEpsilonPEqf_.boundaryField()[patchI][faceI] = 0.0;
                    }
                    else
                    {
                        // yielding

                        // total equivalent strain matrix
                        // (t is start of time-step)
                        const scalar ebart =
                            oldEpsilonPEqf.boundaryField()[patchI][faceI];
                        // yield stress at start of time-step
                        const scalar s00 =
                            sigmaYf_.boundaryField()[patchI][faceI];
                        scalar s0 = 0.0; // updated in loop below
                        const scalar qe =
                            newSigmaEqElasticf.boundaryField()[patchI][faceI];
                        const scalar G = muf.boundaryField()[patchI][faceI];

                        // Calculate deq and s0
                        aravasNewtonLoop
                        (
                            DEpsilonPEqf_.boundaryField()[patchI][faceI],
                            s0,
                            //s00,
                            G,
                            ebart,
                            qe,
                            faceCells[faceI],
                            maxMagDEpsilon
                        );

                        // Update increment of yield stress
                        DSigmaYf_.boundaryField()[patchI][faceI] = s0 - s00;
                    }
                }
            }
        }

        DSigmaYf_.correctBoundaryConditions();
        DEpsilonPEqf_.correctBoundaryConditions();
        plasticNf_.correctBoundaryConditions();

        // Update plastic strain increment
        DEpsilonPf_ = DEpsilonPEqf_*plasticNf_;

        // Relax plastic strain increment to help convergence
        DEpsilonPf_.relax();
    }
}


void newAravasMises::updateYieldStress()
{
    Info << nl << "Updating the yield stress" << endl;
    sigmaY_ += DSigmaY_;
//     sigmaYf_ += DSigmaYf_;

    Info << "\tMax DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;
    epsilonPEq_ += DEpsilonPEq_;

    epsilonP_ += DEpsilonP_;

//     Info << "\tMax DEpsilonPEqf is " << gMax(DEpsilonPEq_) << endl;
//     epsilonPEqf_ += DEpsilonPEqf_;

    // count cells actively yielding
    int numCellsYielding = 0;

    forAll(activeYield_.internalField(), celli)
    {
        if( DEpsilonPEq_.internalField()[celli] > SMALL)
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
        if(!activeYield_.boundaryField()[patchi].coupled())
        {
            forAll(activeYield_.boundaryField()[patchi], facei)
            {
                if( DEpsilonPEq_.boundaryField()[patchi][facei] > SMALL)
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
}

} // end of namespace
// ************************************************************************* //
