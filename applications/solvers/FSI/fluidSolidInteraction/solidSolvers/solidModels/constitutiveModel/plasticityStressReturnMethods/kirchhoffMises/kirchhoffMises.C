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

#include "kirchhoffMises.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "constitutiveModel.H"
#include "transformGeometricField.H"
#include "fvc.H"
#include "solidContactFvPatchVectorField.H"
#include "expElasticPlastic.H"

#include "solidSolver.H"
#include "ULLSMaterialInterface.H"
#include "multiMaterial.H"

#include "processorFvsPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(kirchhoffMises, 0);
    addToRunTimeSelectionTable
    (
        plasticityStressReturn,
        kirchhoffMises,
        dictionary
    );


    // Static variables

        // Tolerance for Newton loop
        scalar kirchhoffMises::LoopTol_ = 1e-8;

        // Maximum number of iterations for Newton loop
        label kirchhoffMises::MaxNewtonIter_ = 100;

        // finiteDiff is the delta for finite difference differentiation
        scalar kirchhoffMises::finiteDiff_ = 0.25e-6;

        // Store sqrt(2/3) as we use it often
        scalar kirchhoffMises::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// yield stress
scalar kirchhoffMises::curYieldStress
(
    const scalar curEpsilonPEq, // current plastic equivalent strain
    const scalar J, // Current Jacobian
    const label cellID // cellID is needed is case in multiMaterial
) const
{
    // We assume that the stress-strain curve was specifed as Cauchy stress vs
    // true strain, but we want the Kirchhoff (tau) yield stress,
    // so we multiply Cauchy stress by J as tauSigmaY = J*sigmaCauchySigmaY

    return J*constitutiveModel_.sigmaY
        (
            max(curEpsilonPEq, SMALL),
            cellID
        );
    // return constitutiveModel_.sigmaY
    //     (
    //         max(curEpsilonPEq, SMALL),
    //         cellID
    //     );
}


scalar kirchhoffMises::yieldFunction
(
    const scalar magSTrial,
    const scalar DLambda,
    const scalar muBar,
    const scalar J,
    const label cellID
) const
{
    // Evaluate current yield function
    // fy = mag(s) - sqrt(2/3)*curSigmaY
    // fy = mag(sTrial - 2*muBar*DLambda*plasticN) - ::sqrt(2.0/3.0)*curSigmaY;
    // fy = magSTrial - 2*muBar*DLambda - ::sqrt(2.0/3.0)*curSigmaY;
    // where
    // fy is the current value of the yield function - zero at convergence.
    // s is the current deviatoric component of tau
    // sTrial is trial version of s
    // plasticN is the return direction
    // DLambda is the current increment of plastic strain multiplier
    // curSigmaY is the current tau yield stress which is typically a function
    // of total equivalent plastic strain (epsilonPEq + DEpsilonPEq)

    return
        magSTrial - 2*muBar*DLambda
      - sqrtTwoOverThree_
       *curYieldStress
        (
            epsilonPEq_[cellID] + sqrtTwoOverThree_*DLambda,
            J,
            cellID
        );
}


void kirchhoffMises::newtonLoop
(
    scalar& DLambda,
    scalar& curSigmaY,
    const scalar magSTrial,
    const scalar muBar,
    const scalar J,
    const label cellID,
    const scalar maxMagDEpsilon
) const
{
    // Loop to determine DEpsilonPEq
    // using Newtion's method

    int i = 0;
    scalar fTrial = yieldFunction(magSTrial, DLambda, muBar, J, cellID);
    scalar residual = 1.0;
    do
    {
        // Numerically calculate derivative of yield function fTrial

        // First Order
        // Hauser 2009 suggested First Order is more efficient for Newton's
        // Method as only two function evaluaitons are required.

        // fTrial step is the the value of fTrial after a small finite
        // difference step
        const scalar fTrialStep  =
            yieldFunction(magSTrial, DLambda + finiteDiff_, muBar,  J, cellID);
        // Numerical derivative of fTrial
        const scalar fTrialDerivative = (fTrialStep - fTrial)/finiteDiff_;

        // Update DLambda
        residual = fTrial/fTrialDerivative;
        DLambda -= residual;

        residual /= maxMagDEpsilon; // Normalise wrt max strain increment

        // fTrial will go to zero at convergence
        fTrial = yieldFunction(magSTrial, DLambda, muBar,  J, cellID);

        if (i == MaxNewtonIter_)
        {
            Warning
                << "Plasticity Newton loop not converging, fx is "
                << mag(fTrial) << endl;
        }
    }
    while ((mag(residual) > LoopTol_) && ++i < MaxNewtonIter_);

    // Update current yield stress
    curSigmaY =
        curYieldStress
        (
            epsilonPEq_[cellID] + sqrtTwoOverThree_*DLambda,
            J,
            cellID
        );
}


void kirchhoffMises::newtonLoop
(
    scalar& DLambda,
    scalar& curSigmaY,
    const scalar magSTrial,
    const scalar muBar,
    const scalar epsilonPEq,
    const scalar J,
    const label cellID
) const
{
    scalar delta = 1;
    label i = 0;
    scalar residual = GREAT;

//     scalar oldSigmaY = curSigmaY;

    // Initialise DLambda
    DLambda = 0;

    scalar curEpsilonPEq;

//     // Reference to exponential hardening law
//     const expElasticPlastic& expLaw =
//         refCast<const expElasticPlastic>(constitutiveModel_.law());

    const rheologyLaw& law = constitutiveModel_.law();

    do
    {
        curEpsilonPEq = epsilonPEq + sqrtTwoOverThree_*DLambda;

        scalar curF =
            (magSTrial - 2*muBar*DLambda)/J
          - sqrtTwoOverThree_*law.sigmaY(curEpsilonPEq, cellID);

        scalar curDF =
           -2*muBar/J - (2.0/3.0)*law.dSigmaY(curEpsilonPEq, cellID);
//            -2*muBar - (2.0/3.0)*J*expLaw.dSigmaY(curEpsilonPEq, 0);

        scalar prevDLambda = DLambda;
        DLambda -= delta*curF/curDF;

        residual = mag(DLambda - prevDLambda);
    }
    while( (residual > LoopTol_) && (++i < MaxNewtonIter_) );

    curEpsilonPEq = epsilonPEq + sqrtTwoOverThree_*DLambda;

    curSigmaY = law.sigmaY(curEpsilonPEq, cellID);

//     curSigmaY = J*expLaw.sigmaY(curEpsilonPEq, 0);

//     scalar DSigmaY = curSigmaY - oldSigmaY;
//     if (DSigmaY < -SMALL)
//     {
//         Info << i << ", " << residual << ", " << DLambda << ", "
//             << DSigmaY << ", "
//             << curSigmaY  << ", " << curEpsilonPEq << ", " << J << endl;
//     }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
kirchhoffMises::kirchhoffMises
(
    const word& name,
    constitutiveModel& constitutiveModel
)
    :
    plasticityStressReturn(name, constitutiveModel),
    constitutiveModel_(constitutiveModel),
    restarted_
    (
        IOobject
        (
            "sigmaY",
            constitutiveModel.sigma().time().timeName(),
            constitutiveModel.sigma().mesh(),
            IOobject::MUST_READ
        ).headerOk()
    ),
    sigmaY_
    (
        IOobject
        (
            "sigmaY",
            constitutiveModel.sigma().time().timeName(),
            constitutiveModel.sigma().mesh(),
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
            constitutiveModel.sigma().mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(sigmaY_, "sigmaY")
    ),
    DSigmaY_
    (
        IOobject
        (
            "DSigmaY",
            sigmaY_.time().timeName(),
            sigmaY_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
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
        dimensionedSymmTensor("zero", dimless, symmTensor::zero) //,
        //zeroGradientFvPatchSymmTensorField::typeName
    ),
    DEpsilonPf_
    (
        IOobject
        (
            "DEpsilonPf",
            sigmaY_.time().timeName(),
            sigmaY_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
            // IOobject::NO_READ,
            // IOobject::NO_WRITE
        ),
        sigmaY_.mesh(),
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
//     bEbar_
//     (
//         IOobject
//         (
//             "bEbar",
//             sigmaY_.time().timeName(),
//             sigmaY_.db(),
//             IOobject::READ_IF_PRESENT,
//             IOobject::AUTO_WRITE
//         ),
//         sigmaY_.mesh(),
//         dimensionedSymmTensor("I", dimless, I)
//     ),
//     bEbarf_
//     (
//         IOobject
//         (
//             "bEbarf",
//             sigmaY_.time().timeName(),
//             sigmaY_.db(),
//             IOobject::READ_IF_PRESENT,
//             IOobject::AUTO_WRITE
//         ),
//         sigmaY_.mesh(),
//         dimensionedSymmTensor("I", dimless, I)
//     ),
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
    DLambda_
    (
        IOobject
        (
            "DLambda",
            sigmaY_.time().timeName(),
            sigmaY_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sigmaY_.mesh(),
        dimensionedScalar("0", dimless, 0.0)
    ),
    DLambdaf_
    (
        IOobject
        (
            "DLambdaf",
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
    ),
    nonLinearPlasticity_(constitutiveModel.law().nonLinearPlasticity()),
    HPtr_(NULL),
    HfPtr_(NULL),
    extrapBoundary_(sigmaY_.boundaryField().size(), false)
{
    Info<< "Creating kirchhoffMises stress return method" << endl;

    if (constitutiveModel_.law().type() == multiMaterial::typeName)
    {
        const unallocLabelList& own = sigmaY_.mesh().owner();

        // Interpolate sigmaY using owner side extrapolation
        forAll(sigmaYf_, faceI)
        {
            sigmaYf_[faceI] = sigmaY_[own[faceI]];
        }

        forAll(sigmaYf_.boundaryField(), patchI)
        {
            sigmaYf_.boundaryField()[patchI] =
                sigmaY_.boundaryField()[patchI].patchInternalField();

//             if
//             (
//                 isA<processorPolyPatch>
//                 (
//                     sigmaY_.mesh().boundaryMesh()[patchI]
//                 )
//             )
//             {
//                 const processorPolyPatch & procPatch =
//                     refCast<const processorPolyPatch>
//                     (
//                         sigmaY_.mesh().boundaryMesh()[patchI]
//                     );

//                 if (procPatch.owner())
//                 {
//                     sigmaYf_.boundaryField()[patchI] =
//                         sigmaY_.boundaryField()[patchI]
//                        .patchInternalField();
//                 }
//                 else
//                 {
//                     sigmaYf_.boundaryField()[patchI] =
//                         sigmaY_.boundaryField()[patchI];
//                 }
//             }
        }
    }

//     // Looking up solid solver
//     const solidSolver& solid =
//         sigmaY_.mesh().lookupObject<solidSolver>
//         (
//             "solidProperties"
//         );

//     if (solid.interface().valid())
//     {
//         sigmaYf_ = solid.interface()->interpolate(sigmaY_);

//         // Initialise interface sigmaY
//         if (isA<ULLSMaterialInterface>(solid.interface()()))
//         {
//             ULLSMaterialInterface& interface =
//                 const_cast<ULLSMaterialInterface&>
//                 (
//                     refCast<const ULLSMaterialInterface>
//                     (
//                         solid.interface()()
//                     )
//                 );

//             const unallocLabelList& ngb = sigmaY_.mesh().neighbour();
//             const unallocLabelList& own = sigmaY_.mesh().owner();

//             forAll(interface.faces(), faceI)
//             {
//                 label curFace = interface.faces()[faceI];

//                 // Internal faces
//                 if (curFace < sigmaY_.mesh().nInternalFaces())
//                 {
//                     label ownCell = own[curFace];
//                     label ngbCell = ngb[curFace];

//                     interface.ownSigmaY()[faceI] =
//                         constitutiveModel_.law().sigmaY(0, ownCell);

//                     interface.ngbSigmaY()[faceI] =
//                         constitutiveModel_.law().sigmaY(0, ngbCell);

//                     Info << interface.ownSigmaY()[faceI] << ", "
//                         << interface.ngbSigmaY()[faceI] << endl;
//                 }
//                 else
//                 {
//                     // Todo
//                 }
//             }
//         }

//     }

    if (nonLinearPlasticity_)
    {
        Info<< "    Plasticity is nonlinear" << endl;
    }
    else
    {
        // Define plastic modulus
        HPtr_ = new volScalarField(constitutiveModel.law().Ep());
        HfPtr_ = new surfaceScalarField(fvc::interpolate(*HPtr_, "Ep"));

        // Check for perfect plasticity
        if (mag(gMax(*HPtr_)) < SMALL)
        {
            Info<< "    Perfect Plasticity" << endl;
            delete HPtr_;
            HPtr_ = NULL;
            delete HfPtr_;
            HfPtr_ = NULL;
        }
        else
        {
            Info<< "    Plasticity is linear" << endl;
        }
    }

    // On restart, the values on symmetryPlanes of surface fields are all zero
    // for some reason. So we will reset them here from the equivalent volume
    // fields
    const fvMesh& mesh = sigmaY_.mesh();
    if (restarted_)
    {
        Info<< "Restarted: plasticity model correcting surface boundary values"
            << endl;
        forAll(mesh.boundary(), patchI)
        {
            if (mesh.boundary()[patchI].type() == "symmetryPlane")
            {
                sigmaYf_.boundaryField()[patchI] =
                    sigmaY_.boundaryField()[patchI];
                DEpsilonPf_.boundaryField()[patchI] =
                    DEpsilonP_.boundaryField()[patchI];
//                 bEbarf_.boundaryField()[patchI] =
//                     bEbar_.boundaryField()[patchI];
                DEpsilonPEqf_.boundaryField()[patchI] =
                    DEpsilonPEq_.boundaryField()[patchI];
                DLambdaf_.boundaryField()[patchI] =
                    DLambda_.boundaryField()[patchI];
                epsilonPEqf_.boundaryField()[patchI] =
                    epsilonPEq_.boundaryField()[patchI];
            }
        }

        // Check if sigmaY is zero anywhere
        // This may happen after mergeMeshesAndFields as the new roller sigmaY
        // field will default to zero
        Info<< "Restarted: checking for zero yield stress"
            << endl;
        scalarField& sigmaYI = sigmaY_.internalField();
        scalarField& sigmaYfI = sigmaYf_.internalField();

        forAll(sigmaYI, cellI)
        {
            if (sigmaYI[cellI] < SMALL)
            {
                sigmaYI[cellI] = GREAT;
            }
        }

        forAll(sigmaYfI, faceI)
        {
            if (sigmaYfI[faceI] < SMALL)
            {
                sigmaYfI[faceI] = GREAT;
            }
        }

        forAll(sigmaY_.boundaryField(), patchI)
        {
            forAll(sigmaY_.boundaryField()[patchI], faceI)
            {
                if (sigmaY_.boundaryField()[patchI][faceI] < SMALL)
                {
                    sigmaY_.boundaryField()[patchI][faceI] = GREAT;
                    sigmaYf_.boundaryField()[patchI][faceI] = GREAT;
                }
            }
        }
        Info<< "    writing sigmaY" << endl;
        sigmaY_.write();
        sigmaYf_.write();
    }

    // Set contact patches to be extrapolated
    // const volVectorField& DU =
    //     mesh.objectRegistry::lookupObject<volVectorField>("DU");
    // forAll(extrapBoundary_, patchI)
    // {
    //     if
    //     (
    //         DU.boundaryField()[patchI].type() ==
    //             solidContactFvPatchVectorField::typeName
    //     )
    //     {
    //         Info<< "    plasticity will be extrapolated on patch "
    //             << mesh.boundary()[patchI].name() << endl;
    //         extrapBoundary_[patchI] = true;
    //     }
    // }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

kirchhoffMises::~kirchhoffMises()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kirchhoffMises::correct()
{
    correctCell();
}

void kirchhoffMises::correctf()
{
    correctFace();
}

void kirchhoffMises::correctCell()
{
    const fvMesh& mesh = sigmaY_.mesh();

    // Elastic properties
    const volScalarField& mu =
        mesh.objectRegistry::lookupObject<volScalarField>("mu");

    // Compute elastic predictor

    // Lookup relative deformation gradient
    // const volTensorField& relF =
    //     mesh.objectRegistry::lookupObject<volTensorField>("relF");
    // volTensorField relFbar = pow(det(relF), -1.0/3.0)*relF;
    // Lookup relative deformation gradient bar
    //const volTensorField& relFbar =
    //  mesh.objectRegistry::lookupObject<volTensorField>("relFbar");
    //volTensorField bEbarTrial =
    //    symm(relFbar & bEbar_.oldTime() & relFbar.T());
    //volSymmTensorField bEbarTrial = transform(relFbar, bEbar_.oldTime());
    const volSymmTensorField& bEbarTrial =
          mesh.objectRegistry::lookupObject<volSymmTensorField>("bBar");

    volSymmTensorField sTrial = mu*dev(bEbarTrial);
    const volScalarField& J =
        mesh.objectRegistry::lookupObject<volScalarField>("J");

    // Calculate trial equivalent stress
    //const volScalarField sigmaEqTrial = sqrt((3.0/2.0)*magSqr(sTrial));

    const volScalarField Ibar = tr(bEbarTrial)/3.0;
    const volScalarField muBar = Ibar*mu;

    // Check for plastic loading
    // and calculate increment of plastic equivalent strain
    // i.e. the plastic multiplier

    // Store previous iteration for under-relaxation and
    // calculation of plastic residual in the solver
    DEpsilonP_.storePrevIter();

    // Normalise residual in Newton method with respect to mag(bE)

//     const scalar maxMagBE = max(gMax(mag(bEbarTrial.internalField())), SMALL);
//     const scalar maxMagBE = max(gMax(mag(bEbar_.internalField())), SMALL);

    // Trial yield function
    // sigmaY is the Cauchy yield stress so we scale it by J
//     const volScalarField fTrial = mag(sTrial) - sqrtTwoOverThree_*J*sigmaY_;
    // sigmaY_ is already scaled by J
    const volScalarField fTrial = mag(sTrial)/J - sqrtTwoOverThree_*sigmaY_;

    // Calculate DLambda_ and plasticN_
    forAll (DEpsilonP_, cellI)
    {
        // Calculate return direction plasticN
        const scalar magS = mag(sTrial[cellI]);
        if (magS > SMALL)
        {
            plasticN_[cellI] = sTrial[cellI]/magS;
        }

        // Calculate DLambda/DEpsilonPEq
        if (fTrial[cellI] < SMALL)
        {
            // elastic
            DSigmaY_[cellI] = 0.0;
            DLambda_[cellI] = 0.0;
        }
        else
        {
            if (nonLinearPlasticity_)
            {
                // Total equivalent plastic strain
                // where t is start of time-step
                scalar curSigmaY = sigmaY_[cellI]; // updated in loop below
//                 scalar curSigmaY = 0.0; // updated in loop below

                // Calculates DEpsilonPEq using Newtons's method
                newtonLoop
                (
                    DLambda_[cellI],
                    curSigmaY,
                    mag(sTrial[cellI]),
                    muBar[cellI],
                    epsilonPEq_[cellI],
                    J[cellI],
                    cellI
                );

//                 if (isA<expElasticPlastic>(constitutiveModel_.law()))
//                 {
//                     // Calculates DEpsilonPEq using Newtons's method
//                     newtonLoop
//                     (
//                         DLambda_[cellI],
//                         curSigmaY,
//                         mag(sTrial[cellI]),
//                         muBar[cellI],
//                         epsilonPEq_[cellI],
//                         J[cellI],
//                         cellI
//                     );
//                 }
//                 else
//                 {
//                     // Calculates DEpsilonPEq using Newtons's method
//                     newtonLoop
//                     (
//                         DLambda_[cellI],
//                         curSigmaY,
//                         mag(sTrial[cellI]),
//                         muBar[cellI],
//                         J[cellI],
//                         cellI,
//                         maxMagBE
//                     );
//                 }

                // Update increment of yield stress
                DSigmaY_[cellI] = curSigmaY - sigmaY_[cellI];
            }
            else
            {
                // Plastic modulus is linear
                DLambda_[cellI] = fTrial[cellI]/(2*muBar[cellI]);

                if (HPtr_)
                {
                    DLambda_[cellI] /=
                        1.0 + (*HPtr_)[cellI]/(3*muBar[cellI]);

                    // Update increment of yield stress
                    DSigmaY_[cellI] = DLambda_[cellI]*(*HPtr_)[cellI];
                }
            }
        }
    }

    forAll (DEpsilonP_.boundaryField(), patchI)
    {
        if
        (
            !DEpsilonP_.boundaryField()[patchI].coupled()
         && !extrapBoundary_[patchI]
        )
        {
            const labelList& faceCells = mesh.boundary()[patchI].faceCells();

            forAll (DEpsilonP_.boundaryField()[patchI], faceI)
            {
                // Calculate direction plasticN
                const scalar magS = mag(sTrial.boundaryField()[patchI][faceI]);
                if (magS > SMALL)
                {
                    plasticN_.boundaryField()[patchI][faceI] =
                        sTrial.boundaryField()[patchI][faceI]/magS;
                }

                // Calculate DEpsilonPEq
                if (fTrial.boundaryField()[patchI][faceI] < SMALL)
                {
                    // elasticity
                    DSigmaY_.boundaryField()[patchI][faceI] = 0.0;
                    DLambda_.boundaryField()[patchI][faceI] = 0.0;
                }
                else
                {
                    // yielding
                    if (nonLinearPlasticity_)
                    {
                        scalar curSigmaY =
                            sigmaY_.boundaryField()[patchI][faceI];
//                         scalar curSigmaY = 0.0; // updated in loop below

                        // Calculate DEpsilonPEq and curSigmaY
                        newtonLoop
                        (
                            DLambda_.boundaryField()[patchI][faceI],
                            curSigmaY,
                            mag(sTrial.boundaryField()[patchI][faceI]),
                            muBar.boundaryField()[patchI][faceI],
                            epsilonPEq_.boundaryField()[patchI][faceI],
                            J.boundaryField()[patchI][faceI],
                            faceCells[faceI]
                        );

//                         if (isA<expElasticPlastic>(constitutiveModel_.law()))
//                         {
//                             // Calculate DEpsilonPEq and curSigmaY
//                             newtonLoop
//                             (
//                                 DLambda_.boundaryField()[patchI][faceI],
//                                 curSigmaY,
//                                 mag(sTrial.boundaryField()[patchI][faceI]),
//                                 muBar.boundaryField()[patchI][faceI],
//                                 epsilonPEq_.boundaryField()[patchI][faceI],
//                                 J.boundaryField()[patchI][faceI],
//                                 faceCells[faceI]
//                             );
//                         }
//                         else
//                         {
//                             // Calculate DEpsilonPEq and curSigmaY
//                             newtonLoop
//                             (
//                                 DLambda_.boundaryField()[patchI][faceI],
//                                 curSigmaY,
//                                 mag(sTrial.boundaryField()[patchI][faceI]),
//                                 muBar.boundaryField()[patchI][faceI],
//                                 J.boundaryField()[patchI][faceI],
//                                 faceCells[faceI],
//                                 maxMagBE
//                             );
//                         }

                        // Update increment of yield stress
                        DSigmaY_.boundaryField()[patchI][faceI] =
                            curSigmaY - sigmaY_.boundaryField()[patchI][faceI];
                    }
                    else
                    {
                        // Plastic modulus is linear
                        DLambda_.boundaryField()[patchI][faceI] =
                            fTrial.boundaryField()[patchI][faceI]
                           /(2*muBar.boundaryField()[patchI][faceI]);

                        if (HPtr_)
                        {
                            DLambda_.boundaryField()[patchI][faceI] /=
                                1.0 + (*HPtr_).boundaryField()[patchI][faceI]
                               /(3*muBar.boundaryField()[patchI][faceI]);

                            // Update increment of yield stress
                            DSigmaY_.boundaryField()[patchI][faceI] =
                                DLambda_.boundaryField()[patchI][faceI]
                                *(*HPtr_).boundaryField()[patchI][faceI];
                        }
                    }
                }
            }
        }
    }

    DSigmaY_.correctBoundaryConditions();
    DLambda_.correctBoundaryConditions();
    plasticN_.correctBoundaryConditions();

    // Extrapolating from the internal field on boundaries that cause
    // convergence problems
    forAll(DEpsilonP_.boundaryField(), patchI)
    {
        if (extrapBoundary_[patchI])
        {
            DSigmaY_.boundaryField()[patchI] =
                DSigmaY_.boundaryField()[patchI].patchInternalField();

            DLambda_.boundaryField()[patchI] =
                DLambda_.boundaryField()[patchI].patchInternalField();

            plasticN_.boundaryField()[patchI] =
                plasticN_.boundaryField()[patchI].patchInternalField();
        }
    }

    // Update DEpsilonP and DEpsilonPEq
    DEpsilonPEq_ = sqrtTwoOverThree_*DLambda_;
    DEpsilonP_ = Ibar*DLambda_*plasticN_;
    DEpsilonP_.relax();

//     if (true)
//     {
//         this->correctf();
//     }

//     // Update bEbar
//     volSymmTensorField s = sTrial - 2*mu*DEpsilonP_;
//     bEbar_ = (s/mu) + Ibar*I;
}


void kirchhoffMises::correctFace()
{
    // This is the same procedure as in correct() except DEpsilonPf is
    // calculated instead of DEpsilonP
    const fvMesh& mesh = sigmaYf_.mesh();

    // Elastic properties
    const surfaceScalarField& muf =
        mesh.objectRegistry::lookupObject<surfaceScalarField>("muf");

    // Compute elastic predictor
    const surfaceSymmTensorField& bEbarTrialf =
          mesh.objectRegistry::lookupObject<surfaceSymmTensorField>("bBarf");

    const surfaceSymmTensorField sTrialf = muf*dev(bEbarTrialf);
    const surfaceScalarField magSTrialf = mag(sTrialf);
    const surfaceScalarField& Jf =
        mesh.objectRegistry::lookupObject<surfaceScalarField>("Jf");

    const surfaceScalarField Ibarf = tr(bEbarTrialf)/3.0;
    const surfaceScalarField muBarf = Ibarf*muf;

    // Check for plastic loading - plane direction and
    // plastic scalar multiplier

    // Store previous iteration for under-relaxation and calculation of plastic
    // residual in the solver
    DEpsilonPf_.storePrevIter();

    // Normalise residual in Newton method with respect to mag(bE)
//     const scalar maxMagBE =
//         max(gMax(mag(bEbarTrialf.internalField())), SMALL);

//     const scalar maxMagBE = max(gMax(mag(bEbarf_.internalField())), SMALL);

    // Trial yield function
    // sigmaY is the Cauchy yield stress so we scale it by J

//     const surfaceScalarField fTrialf =
//         mag(sTrialf) - sqrtTwoOverThree_*Jf*sigmaYf_;

    const surfaceScalarField fTrialf =
        mag(sTrialf)/Jf - sqrtTwoOverThree_*sigmaYf_;

    const unallocLabelList& own = mesh.owner();

    // Calculate DLambdaf_ and plasticNf_
    forAll (DEpsilonPf_, faceI)
    {
        // Calculate return direction plasticN
        if (magSTrialf[faceI] > SMALL)
        {
            plasticNf_[faceI] = sTrialf[faceI]/magSTrialf[faceI];
        }

        // Calculate DLambda/DEpsilonPEq
        if (fTrialf[faceI] < SMALL)
        {
            // elastic
            DSigmaYf_[faceI] = 0.0;
            DLambdaf_[faceI] = 0.0;
        }
        else
        {
            if (nonLinearPlasticity_)
            {
                // Total equivalent plastic strain where t is start
                // of time-step
                scalar curSigmaYf = sigmaYf_[faceI];
//                 scalar curSigmaYf = 0.0; // updated in loop below

                // Calculates DEpsilonPEq using Newtons's method
                newtonLoop
                (
                    DLambdaf_[faceI],
                    curSigmaYf,
                    magSTrialf[faceI],
                    muBarf[faceI],
                    epsilonPEqf_[faceI],
                    Jf[faceI],
                    own[faceI]
                );

//                 if (isA<expElasticPlastic>(constitutiveModel_.law()))
//                 {
//                     // Calculates DEpsilonPEq using Newtons's method
//                     newtonLoop
//                     (
//                         DLambdaf_[faceI],
//                         curSigmaYf,
//                         magSTrialf[faceI],
//                         muBarf[faceI],
//                         epsilonPEqf_[faceI],
//                         Jf[faceI],
//                         own[faceI]
//                     );
//                 }
//                 else
//                 {
//                     // Calculates DEpsilonPEq using Newtons's method
//                     newtonLoop
//                     (
//                         DLambdaf_[faceI],
//                         curSigmaYf,
//                         magSTrialf[faceI],
//                         muBarf[faceI],
//                         Jf[faceI],
//                         own[faceI],
//                         maxMagBE
//                     );
//                 }

                // Update increment of yield stress
                DSigmaYf_[faceI] = curSigmaYf - sigmaYf_[faceI];
            }
            else
            {
                // Plastic modulus is linear
                DLambdaf_[faceI] = fTrialf[faceI]/(2*muBarf[faceI]);

                if (HfPtr_)
                {
                    DLambdaf_[faceI] /=
                        1.0 + (*HfPtr_)[faceI]/(3*muBarf[faceI]);

                    // Update increment of yield stress
                    DSigmaYf_[faceI] = DLambdaf_[faceI]*(*HfPtr_)[faceI];
                }
            }
        }
    }

    forAll (DEpsilonPf_.boundaryField(), patchI)
    {
        //if (!DEpsilonPf_.boundaryField()[patchI].coupled()) // CHECK
        {
            const labelList& faceCells = mesh.boundary()[patchI].faceCells();

            forAll (DEpsilonPf_.boundaryField()[patchI], faceI)
            {
                // Calculate direction plasticN
                if (magSTrialf.boundaryField()[patchI][faceI] > SMALL)
                {
                    plasticNf_.boundaryField()[patchI][faceI] =
                        sTrialf.boundaryField()[patchI][faceI]
                       /magSTrialf.boundaryField()[patchI][faceI];
                }

                // Calculate DEpsilonPEq
                if (fTrialf.boundaryField()[patchI][faceI] < SMALL)
                {
                    // elasticity
                    DSigmaYf_.boundaryField()[patchI][faceI] = 0.0;
                    DLambdaf_.boundaryField()[patchI][faceI] = 0.0;
                }
                else
                {
                    // yielding
                    if (nonLinearPlasticity_)
                    {
                        scalar curSigmaYf =
                            sigmaYf_.boundaryField()[patchI][faceI];

//                         scalar curSigmaYf = 0.0; // updated in loop below

                        // Calculate DEpsilonPEq and curSigmaY
                        newtonLoop
                        (
                            DLambdaf_.boundaryField()[patchI][faceI],
                            curSigmaYf,
                            mag(sTrialf.boundaryField()[patchI][faceI]),
                            muBarf.boundaryField()[patchI][faceI],
                            epsilonPEqf_.boundaryField()[patchI][faceI],
                            Jf.boundaryField()[patchI][faceI],
                            faceCells[faceI]
                        );

//                         if (isA<expElasticPlastic>(constitutiveModel_.law()))
//                         {
//                             // Calculate DEpsilonPEq and curSigmaY
//                             newtonLoop
//                             (
//                                 DLambdaf_.boundaryField()[patchI][faceI],
//                                 curSigmaYf,
//                                 mag(sTrialf.boundaryField()[patchI][faceI]),
//                                 muBarf.boundaryField()[patchI][faceI],
//                                 epsilonPEqf_.boundaryField()[patchI][faceI],
//                                 Jf.boundaryField()[patchI][faceI],
//                                 faceCells[faceI]
//                             );
//                         }
//                         else
//                         {
//                             // Calculate DEpsilonPEq and curSigmaY
//                             newtonLoop
//                             (
//                                 DLambdaf_.boundaryField()[patchI][faceI],
//                                 curSigmaYf,
//                                 mag(sTrialf.boundaryField()[patchI][faceI]),
//                                 muBarf.boundaryField()[patchI][faceI],
//                                 Jf.boundaryField()[patchI][faceI],
//                                 faceCells[faceI],
//                                 maxMagBE
//                             );
//                         }

                        // Update increment of yield stress
                        DSigmaYf_.boundaryField()[patchI][faceI] =
                            curSigmaYf
                          - sigmaYf_.boundaryField()[patchI][faceI];
                    }
                    else
                    {
                        // Plastic modulus is linear
                        DLambdaf_.boundaryField()[patchI][faceI] =
                            fTrialf.boundaryField()[patchI][faceI]
                           /(2*muBarf.boundaryField()[patchI][faceI]);

                        if (HfPtr_)
                        {
                            DLambdaf_.boundaryField()[patchI][faceI] /=
                                1.0 + (*HfPtr_).boundaryField()[patchI][faceI]
                               /(3*muBarf.boundaryField()[patchI][faceI]);

                            // Update increment of yield stress
                            DSigmaYf_.boundaryField()[patchI][faceI] =
                                DLambdaf_.boundaryField()[patchI][faceI]
                               *(*HfPtr_).boundaryField()[patchI][faceI];
                        }
                    }
                }
            }
        }
    }

    DSigmaYf_.correctBoundaryConditions();
    DLambdaf_.correctBoundaryConditions();
    plasticNf_.correctBoundaryConditions();

    // Update DEpsilonP and DEpsilonPEq
    DEpsilonPEqf_ = sqrtTwoOverThree_*DLambdaf_;
    DEpsilonPf_ = Ibarf*DLambdaf_*plasticNf_;
    DEpsilonPf_.relax();

    // Looking up solid solver
    const solidSolver& solid =
        mesh.lookupObject<solidSolver>
        (
            "solidProperties"
        );

    // Calculate DLambda and plasticN
    // at the neighbour side of the interface
    // Owner side is already calculated
    if (solid.interface().valid())
    {
        if (isA<ULLSMaterialInterface>(solid.interface()()))
        {
            ULLSMaterialInterface& interface =
                const_cast<ULLSMaterialInterface&>
                (
                    refCast<const ULLSMaterialInterface>
                    (
                        solid.interface()()
                    )
                );

            const unallocLabelList& ngb = mesh.neighbour();
//             const unallocLabelList& own = mesh.owner();

            const volTensorField& gradDD =
                mesh.lookupObject<volTensorField>
                (
                    "grad(" + interface.DD().name() + ')'
                );

            const surfaceTensorField& gradDDf =
                mesh.lookupObject<surfaceTensorField>
                (
                    "grad" + interface.DD().name() + 'f'
                );

            const volScalarField& mu =
                mesh.lookupObject<volScalarField>("mu");

            forAll(interface.faces(), faceI)
            {
                label curFace = interface.faces()[faceI];

                // Internal faces
                if (curFace < mesh.nInternalFaces())
                {
//                     label ownCell = own[curFace];
                    label ngbCell = ngb[curFace];

                    scalarField J =
                        interface.J(faceI, gradDD, gradDDf);
                    scalar ngbJ = J[1];

                    symmTensorField BbarTrial =
                        interface.trialBbar(faceI, gradDD, gradDDf);
                    symmTensor ngbBbarTrial = BbarTrial[1];

                    scalar ngbMu = mu[ngbCell];

                    symmTensor ngbSTrial = ngbMu*dev(ngbBbarTrial);

                    scalar ngbMagSTrial = mag(ngbSTrial);

                    scalar ngbFTrial =
                        ngbMagSTrial/ngbJ
                      - sqrtTwoOverThree_*interface.ngbSigmaY()[faceI];

//                     Info << ngbMagSTrial << endl;
//                     Info << ngbJ << endl;
//                     Info << ngbMagSTrial/ngbJ << endl;
//                     Info << interface.ngbSigmaY()[faceI]
//                         << endl;
//                     Info << ngbFTrial << endl;

                    scalar ngbIbar = tr(ngbBbarTrial)/3.0;
                    scalar ngbMuBar = ngbIbar*ngbMu;

                    symmTensor ngbPlasticN = symmTensor::zero;
                    scalar ngbDLambda = 0;

                    // Calculate return direction plasticN
                    if (ngbMagSTrial > SMALL)
                    {
                        ngbPlasticN = ngbSTrial/ngbMagSTrial;
                    }

                    // Calculate DLambda/DEpsilonPEq
                    if (ngbFTrial < SMALL)
                    {
                        // elastic
                        interface.ngbDSigmaY()[faceI] = 0;
                        ngbDLambda = 0;
                    }
                    else
                    {
                        if (nonLinearPlasticity_)
                        {
                            scalar curNgbSigmaY =
                                interface.ngbSigmaY()[faceI];

                            // Calculates DEpsilonPEq using Newtons's method
                            newtonLoop
                            (
                                ngbDLambda,
                                curNgbSigmaY,
                                ngbMagSTrial,
                                ngbMuBar,
                                interface.ngbSigmaY()[faceI],
                                ngbJ,
                                ngbCell
                            );

                            // Update increment of yield stress
                            interface.ngbDSigmaY()[faceI] =
                                curNgbSigmaY - interface.ngbSigmaY()[faceI];
                        }
                        else
                        {
                            // Plastic modulus is linear
                            ngbDLambda = ngbFTrial/(2*ngbMuBar);

                            if (HPtr_)
                            {
                                ngbDLambda =
                                    1.0 + (*HPtr_)[ngbCell]/(3*ngbMuBar);

                                // Update increment of yield stress
                                interface.ngbDSigmaY()[faceI] =
                                    ngbDLambda*(*HPtr_)[ngbCell];
                            }
                        }
                    }

                    interface.ngbDEpsilonPEq()[faceI] =
                        sqrtTwoOverThree_*ngbDLambda;

                    symmTensor prevNgbDEpsilonP =
                        interface.ngbDEpsilonP()[faceI];

                    interface.ngbDEpsilonP()[faceI] =
                        ngbIbar*ngbDLambda*ngbPlasticN;

                    // Relax DEpsilonP
                    {
                        scalar alpha = 0;

                        if (mesh.solutionDict().relaxField(DEpsilonPf_.name()))
                        {
                            alpha =
                                mesh.solutionDict().fieldRelaxationFactor
                                (
                                    DEpsilonPf_.name()
                                );
                        }

                        if (alpha > 0)
                        {
                            interface.ngbDEpsilonP()[faceI] =
                                prevNgbDEpsilonP
                              + alpha
                               *(
                                   interface.ngbDEpsilonP()[faceI]
                                 - prevNgbDEpsilonP
                                );
                        }
                    }

                    // Owner side is already calculated
                    interface.ownDEpsilonPEq()[faceI] =
                        DEpsilonPEqf_[curFace];
                    interface.ownDEpsilonP()[faceI] =
                        DEpsilonPf_[curFace];
                    interface.ownDSigmaY()[faceI] =
                        DSigmaYf_[curFace];
                }
                else
                {
                    // Processor patch face

                    label curPatch =
                        mesh.boundaryMesh().whichPatch(curFace);
                    label curPatchFace =
                        curFace - mesh.boundaryMesh()[curPatch].start();

                    const processorPolyPatch & procPatch =
                        refCast<const processorPolyPatch>
                        (
                            mesh.boundaryMesh()[curPatch]
                        );

                    const unallocLabelList& faceCells =
                        mesh.boundary()[curPatch].faceCells();

                    if (procPatch.owner())
                    {
                        // Owner side is already calculated
                        interface.ownDEpsilonPEq()[faceI] =
                            DEpsilonPEqf_
                           .boundaryField()[curPatch][curPatchFace];
                        interface.ownDEpsilonP()[faceI] =
                            DEpsilonPf_
                           .boundaryField()[curPatch][curPatchFace];
                        interface.ownDSigmaY()[faceI] =
                            DSigmaYf_.boundaryField()[curPatch][curPatchFace];

                        // Neigbour fields must be calculated
                        // at the neighbour side
                    }
                    else
                    {
                        label ownCell = faceCells[curPatchFace];

                        scalarField J =
                            interface.J(faceI, gradDD, gradDDf);
                        scalar ownJ = J[0];

                        symmTensorField BbarTrial =
                            interface.trialBbar(faceI, gradDD, gradDDf);
                        symmTensor ownBbarTrial = BbarTrial[0];

                        scalar ownMu = mu[ownCell];

                        symmTensor ownSTrial = ownMu*dev(ownBbarTrial);

                        scalar ownMagSTrial = mag(ownSTrial);

                        scalar ownFTrial =
                            ownMagSTrial/ownJ
                          - sqrtTwoOverThree_*interface.ownSigmaY()[faceI];

                        scalar ownIbar = tr(ownBbarTrial)/3.0;
                        scalar ownMuBar = ownIbar*ownMu;

                        symmTensor ownPlasticN = symmTensor::zero;
                        scalar ownDLambda = 0;

                        // Calculate return direction plasticN
                        if (ownMagSTrial > SMALL)
                        {
                            ownPlasticN = ownSTrial/ownMagSTrial;
                        }

                        // Calculate DLambda/DEpsilonPEq
                        if (ownFTrial < SMALL)
                        {
                            // elastic
                            interface.ownDSigmaY()[faceI] = 0;
                            ownDLambda = 0;
                        }
                        else
                        {
                            if (nonLinearPlasticity_)
                            {
                                scalar curOwnSigmaY =
                                    interface.ownSigmaY()[faceI];

                                // Calculates DEpsilonPEq using
                                // Newtons's method
                                newtonLoop
                                (
                                    ownDLambda,
                                    curOwnSigmaY,
                                    ownMagSTrial,
                                    ownMuBar,
                                    interface.ownSigmaY()[faceI],
                                    ownJ,
                                    ownCell
                                );

                                // Update increment of yield stress
                                interface.ownDSigmaY()[faceI] =
                                    curOwnSigmaY
                                  - interface.ownSigmaY()[faceI];
                            }
                            else
                            {
                                // Plastic modulus is linear
                                ownDLambda = ownFTrial/(2*ownMuBar);

                                if (HPtr_)
                                {
                                    ownDLambda =
                                        1.0 + (*HPtr_)[ownCell]/(3*ownMuBar);

                                    // Update increment of yield stress
                                    interface.ownDSigmaY()[faceI] =
                                        ownDLambda*(*HPtr_)[ownCell];
                                }
                            }
                        }

                        interface.ownDEpsilonPEq()[faceI] =
                            sqrtTwoOverThree_*ownDLambda;

                        symmTensor prevOwnDEpsilonP =
                            interface.ownDEpsilonP()[faceI];

                        interface.ownDEpsilonP()[faceI] =
                            ownIbar*ownDLambda*ownPlasticN;

                        // Relax DEpsilonP
                        {
                            scalar alpha = 0;

                            if
                            (
                                mesh.solutionDict().relaxField
                                (
                                    DEpsilonPf_.name()
                                )
                            )
                            {
                                alpha =
                                    mesh.solutionDict().fieldRelaxationFactor
                                    (
                                        DEpsilonPf_.name()
                                    );
                            }

                            if (alpha > 0)
                            {
                                interface.ownDEpsilonP()[faceI] =
                                    prevOwnDEpsilonP
                                  + alpha
                                   *(
                                        interface.ownDEpsilonP()[faceI]
                                      - prevOwnDEpsilonP
                                    );
                            }
                        }

                        DEpsilonPEqf_.boundaryField()[curPatch][curPatchFace] =
                            interface.ownDEpsilonPEq()[faceI];

                        DEpsilonPf_.boundaryField()[curPatch][curPatchFace] =
                            interface.ownDEpsilonP()[faceI];

                        DSigmaYf_.boundaryField()[curPatch][curPatchFace] =
                            interface.ownDSigmaY()[faceI];


//                         // Neighbour side is already calculated
//                         interface.ngbDEpsilonPEq()[faceI] =
//                             DEpsilonPEqf_
//                            .boundaryField()[curPatch][curPatchFace];
//                         interface.ngbDEpsilonP() =
//                             DEpsilonPf_
//                            .boundaryField()[curPatch][curPatchFace];
//                         interface.ngbDSigmaY()[faceI] =
//                             DSigmaYf_.boundaryField()[curPatch][curPatchFace];
                    }

//                     interface.ownDEpsilonPEq()[faceI] =
//                         DEpsilonPEqf_.boundaryField()[curPatch][curPatchFace];
//                     interface.ownDEpsilonP() =
//                         DEpsilonPf_.boundaryField()[curPatch][curPatchFace];
//                     interface.ownDSigmaY()[faceI] =
//                         DSigmaYf_.boundaryField()[curPatch][curPatchFace];

//                     interface.ngbDEpsilonPEq()[faceI] =
//                         ngbProcDEpsilonPEq[curPatch][curPatchFace];
//                     interface.ngbDEpsilonP() =
//                         ngbProcDEpsilonP[curPatch][curPatchFace];
//                     interface.ngbDSigmaY()[faceI] =
//                         ngbProcDSigmaY[curPatch][curPatchFace];
                }
            }

            // Processor neighbour fields

            FieldField<Field, scalar> ngbProcDEpsilonPEq
            (
                mesh.boundary().size()
            );
            FieldField<Field, symmTensor> ngbProcDEpsilonP
            (
                mesh.boundary().size()
            );
            FieldField<Field, scalar> ngbProcDSigmaY
            (
                mesh.boundary().size()
            );

#           include "sendReceiveProcFields.H"

            forAll(interface.faces(), faceI)
            {
                label curFace = interface.faces()[faceI];

                // Internal faces
                if (curFace >= mesh.nInternalFaces())
                {
                    // Processor patch face
                    label curPatch =
                        mesh.boundaryMesh().whichPatch(curFace);
                    label curPatchFace =
                        curFace - mesh.boundaryMesh()[curPatch].start();

                    const processorPolyPatch & procPatch =
                        refCast<const processorPolyPatch>
                        (
                            mesh.boundaryMesh()[curPatch]
                        );

//                     const unallocLabelList& faceCells =
//                         mesh.boundary()[curPatch].faceCells();

                    interface.ngbDEpsilonPEq()[faceI] =
                        ngbProcDEpsilonPEq[curPatch][curPatchFace];

                    interface.ngbDEpsilonP()[faceI] =
                        ngbProcDEpsilonP[curPatch][curPatchFace];

                    interface.ngbDSigmaY()[faceI] =
                        ngbProcDSigmaY[curPatch][curPatchFace];

                    if (!procPatch.owner())
                    {
                        DEpsilonPEqf_.boundaryField()[curPatch][curPatchFace] =
                            interface.ngbDEpsilonPEq()[faceI];

                        DEpsilonPf_.boundaryField()[curPatch][curPatchFace] =
                            interface.ngbDEpsilonP()[faceI];

                        DSigmaYf_.boundaryField()[curPatch][curPatchFace] =
                            interface.ngbDSigmaY()[faceI];

                    }
                }
            }
        }
    }


//     // Update bEbar
//     surfaceSymmTensorField sf = sTrialf - 2*muf*DEpsilonPf_;
//     bEbarf_ = (sf/muf) + Ibarf*I;
}


void kirchhoffMises::updateYieldStress()
{
    correctCell();

    Info << nl << "Updating the yield stress" << endl;
    sigmaY_ += DSigmaY_;
    sigmaYf_ += DSigmaYf_;

//     Info << "maxSigmaY: " << max(sigmaY_) << endl;
//     Info << "minSigmaY: " << min(sigmaY_) << endl;
//     Info << "maxSigmaYf: " << max(sigmaYf_) << endl;
//     Info << "minSigmaYf: " << min(sigmaYf_) << endl;

    Info << "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;
    epsilonPEq_ += DEpsilonPEq_;

    Info << "    Max DEpsilonPEqf is " << gMax(DEpsilonPEqf_) << endl;
    epsilonPEqf_ += DEpsilonPEqf_;


//     Info << "    Max DSigmaYf is " << gMax(DSigmaYf_) << endl;
//     Info << "    Min DSigmaYf is " << gMin(DSigmaYf_) << endl;

//     Info << "    Max sigmaYf is " << gMax(sigmaYf_) << endl;
//     Info << "    Min sigmaYf is " << gMin(sigmaYf_) << endl;

//     Info << "    Max DSigmaY is " << gMax(DSigmaY_) << endl;
//     Info << "    Min DSigmaY is " << gMin(DSigmaY_) << endl;
//     Info << "    Max sigmaY is " << gMax(sigmaY_) << endl;
//     Info << "    Min sigmaY is " << gMin(sigmaY_) << endl;

//     Info << "maxEpsilonPEq: " << max(epsilonPEq_) << endl;
//     Info << "minEpsilonPEq: " << min(epsilonPEq_) << endl;
//     Info << "maxEpsilonPEqf: " << max(epsilonPEqf_) << endl;
//     Info << "minEpsilonPEqf: " << min(epsilonPEqf_) << endl;


    // Set activeYield_ field using DEpsilonPEqf_ field
    const fvMesh& mesh = sigmaY_.mesh();
    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();
    activeYield_ *= 0;

    // Looking up solid solver
    const solidSolver& solid =
        sigmaY_.mesh().lookupObject<solidSolver>
        (
            "solidProperties"
        );

    if (solid.interface().valid())
    {
        ULLSMaterialInterface& interface =
            const_cast<ULLSMaterialInterface&>
            (
                refCast<const ULLSMaterialInterface>
                (
                    solid.interface()()
                )
            );


        forAll(DEpsilonPEqf_, faceI)
        {
            label intFace = findIndex(interface.faces(), faceI);

            if (intFace != -1)
            {
                if (interface.ownDEpsilonPEq()[intFace] > SMALL)
                {
                    activeYield_.internalField()[own[faceI]] = 1.0;
                }

                if (interface.ngbDEpsilonPEq()[intFace] > SMALL)
                {
                    activeYield_.internalField()[nei[faceI]] = 1.0;
                }
            }
            else
            {
                if ( DEpsilonPEqf_.internalField()[faceI] > SMALL)
                {
                    activeYield_.internalField()[own[faceI]] = 1.0;
                    activeYield_.internalField()[nei[faceI]] = 1.0;
                }
            }
        }

        // Todo
        forAll(DEpsilonPEqf_.boundaryField(), patchI)
        {
            const labelList& faceCells = mesh.boundary()[patchI].faceCells();

            forAll(DEpsilonPEqf_.boundaryField()[patchI], faceI)
            {
                label start = mesh.boundaryMesh()[patchI].start();
                label globalFaceI = start + faceI;

                label intFace = findIndex(interface.faces(), globalFaceI);

                if (intFace != -1)
                {
                    if (interface.ownDEpsilonPEq()[intFace] > SMALL)
                    {
                        activeYield_.internalField()[faceCells[faceI]] = 1.0;
                    }
                }
                else
                {
                    if ( DEpsilonPEqf_.boundaryField()[patchI][faceI] > SMALL)
                    {
                        activeYield_.internalField()[faceCells[faceI]] = 1.0;
                        activeYield_.boundaryField()[patchI][faceI] = 1.0;
                    }
                }
            }
        }
    }
    else
    {
        forAll(DEpsilonPEqf_, faceI)
        {
            if ( DEpsilonPEqf_.internalField()[faceI] > SMALL)
            {
                activeYield_.internalField()[own[faceI]] = 1.0;
                activeYield_.internalField()[nei[faceI]] = 1.0;
            }
        }
        forAll(DEpsilonPEqf_.boundaryField(), patchI)
        {
            const labelList& faceCells = mesh.boundary()[patchI].faceCells();
            forAll(DEpsilonPEqf_.boundaryField()[patchI], faceI)
            {
                if ( DEpsilonPEqf_.boundaryField()[patchI][faceI] > SMALL)
                {
                    activeYield_.internalField()[faceCells[faceI]] = 1.0;
                    activeYield_.boundaryField()[patchI][faceI] = 1.0;
                }
            }
        }
    }

    activeYield_.correctBoundaryConditions();

    // count cells actively yielding
    int numCellsYielding = 0;

    //
    activeYield_ *= 0;
    forAll(activeYield_.internalField(), celli)
    {
        if ( DEpsilonPEq_.internalField()[celli] > SMALL)
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
                if ( DEpsilonPEq_.boundaryField()[patchi][facei] > SMALL)
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
    activeYield_.correctBoundaryConditions();


//     forAll(activeYield_.internalField(), celli)
//     {
//         if (activeYield_.internalField()[celli] > SMALL)
//         {
//             numCellsYielding++;
//         }
//     }
//     reduce(numCellsYielding, sumOp<int>());


    Info<< "    " << numCellsYielding << " cells are actively yielding"
        << nl << endl;
}


const volSymmTensorField& kirchhoffMises::DEpsilonP() const
{
    return DEpsilonP_;
}

const surfaceSymmTensorField& kirchhoffMises::DEpsilonPf() const
{
    return DEpsilonPf_;
}

} // end of namespace
// ************************************************************************* //
