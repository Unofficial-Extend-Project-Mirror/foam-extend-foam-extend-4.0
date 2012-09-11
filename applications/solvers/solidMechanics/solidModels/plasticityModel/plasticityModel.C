/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
7
\*---------------------------------------------------------------------------*/

#include "plasticityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(plasticityModel, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasticityModel::plasticityModel
(
    const volTensorField& gradDU,
    const volSymmTensorField& epsilon,
    const volSymmTensorField& sigma
)
:
    rheologyModel(sigma),
    gradDU_(gradDU),
    epsilon_(epsilon),
    plasticityModelCoeffs_(subDict(type() + "Coeffs")),

    beta_
    (
        IOobject
        (
            "beta",
            sigma.time().timeName(),
            sigma.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sigma.mesh(),
        dimensionedScalar("0", dimless, 0)
	),
    sigmaY_
    (
        IOobject
        (
            "sigmaY",
            sigma.time().timeName(),
            sigma.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sigmaY()
    ),
    DEpsilonP_
    (
        IOobject
        (
            "DepsilonP",
            sigma.time().timeName(),
            sigma.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sigma.mesh(),
        dimensionedSymmTensor("0", dimless, symmTensor::zero)
    ),
    mu_
    (
        IOobject
        (
            "mu",
            sigma.time().timeName(),
            sigma.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mu()
    ),
    lambda_
    (
        IOobject
        (
            "lambda",
            sigma.time().timeName(),
            sigma.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        lambda()
    )
{}


plasticityModel::~plasticityModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void plasticityModel::correct()
{
//    rheologyModel::correct();
    Info << "\tCorrecting plasticity model ... " << flush;

    const volSymmTensorField DEpsilon = 
        symm(gradDU_)
      + dimensioned<symmTensor>
        (
            "SMALL",
            dimless,
            symmTensor(SMALL, SMALL, SMALL, SMALL, SMALL, SMALL)
        );

    const volScalarField epsilonEq =
        sqrt((2.0/3.0)*magSqr(dev(epsilon_ + DEpsilon))) 
      + dimensionedScalar("SMALL", dimless, SMALL);

    scalarField& sigmaYI = sigmaY_.internalField();

    scalarField initialSigmaYI = sigmaY()().internalField();

    volScalarField epsilonEqCorr = epsilonEq;
    /*
    forAll(mu_.internalField(), cellI)
    {
        if (sigmaYI[cellI] > initialSigmaYI[cellI])
        {
            epsilonEqCorr.internalField()[cellI] = 0.02;
        }
    }
    */
    // Update mu and lambda
    mu_= mu(epsilonEqCorr);
    lambda_ = lambda(epsilonEqCorr);

    // somewhat underestimating - should be combination/line search of epsEq_old and epsEq!!!
    //    volScalarField Ep_ = Ep(epsilonEq);

    const volScalarField DEpsilonEq =
        sqrt((2.0/3.0)*magSqr(dev(epsilon_ + DEpsilon)))
      - sqrt((2.0/3.0)*magSqr(dev(epsilon_)))
      + dimensionedScalar("SMALL", dimless, SMALL);

    const volSymmTensorField DSigma = 
        2*mu_*(DEpsilon - DEpsilonP_) + I*(lambda_*tr(DEpsilon));

    const volSymmTensorField& oldSigma = sigma();

    const volScalarField oldSigmaEq = sqrt(1.5*magSqr(dev(oldSigma)));

    const volSymmTensorField sigma_ = sigma() + DSigma;

    const volScalarField sigmaEq = 
        sqrt(1.5*magSqr(dev(sigma_)))
      + dimensionedScalar("SMALL", dimPressure, SMALL);

    const volSymmTensorField devSigma = dev(sigma_);

    const volSymmTensorField DSigmaE = DSigma + 2*mu_*DEpsilonP_;

    const volScalarField sigmaEqE  = sqrt(1.5*magSqr(dev(oldSigma + DSigmaE)));

    const volScalarField DSigmaEqE = sqrt(1.5*magSqr(dev(DSigmaE)));


    volScalarField Ep_ = Ep(sigmaEq);

    // Update internal beta
    const scalarField& muI = mu_.internalField();

    const scalarField& lambdaI = lambda_.internalField();

    const symmTensorField& DEpsilonI = DEpsilon.internalField();

    const scalarField& DEpsilonEqI  = DEpsilonEq.internalField();

    const symmTensorField& oldSigmaI = oldSigma.internalField();

    const scalarField& oldSigmaEqI = oldSigmaEq.internalField();

    const symmTensorField& devSigmaI = devSigma.internalField();

    const symmTensorField& DSigmaEI = DSigmaE.internalField();

    const scalarField& sigmaEqEI  = sigmaEqE;

    const scalarField& DSigmaEqEI = DSigmaEqE;

    const scalarField& oldBetaI = beta_.oldTime().internalField();

    scalarField& betaI = beta_.internalField();

    forAll (betaI, cellI)
    {
        tensor curDEpsEPred = tensor::zero;

        if( (DEpsilonEqI[cellI] >= 0) && (oldBetaI[cellI] > SMALL) )
        {
            betaI[cellI] = 1.0;
            curDEpsEPred = tensor::zero;
        }
        else
        {
            betaI[cellI] = 0.0;
            curDEpsEPred = DEpsilonI[cellI];
 
            if
            ( 
                (DEpsilonEqI[cellI] >= 0)
             && (sigmaEqEI[cellI] >= sigmaYI[cellI])
            )
            {
                scalar C = sqr(oldSigmaEqI[cellI]) - sqr(sigmaYI[cellI]);
                scalar B = 3.0*(dev(oldSigmaI[cellI]) && dev(DSigmaEI[cellI]));
                scalar A = sqr(DSigmaEqEI[cellI]);
 
		scalar alpha = (-B + ::sqrt(mag(B*B - 4*A*C)))/(2*A + SMALL);
                //   scalar alpha = (-B + ::sqrt((B*B - 4*A*C)))/(2*A + SMALL); 
                curDEpsEPred =
                    alpha/(2.0*muI[cellI] + SMALL)
                   *(
                        DSigmaEI[cellI] 
                      - (lambdaI[cellI]/(2*muI[cellI] + 3*lambdaI[cellI] + SMALL))
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
        const scalarField& muPatch = mu_.boundaryField()[patchI];
        const scalarField& lambdaPatch = lambda_.boundaryField()[patchI];

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
                betaPatch[faceI] = 1;
                curDEpsEPred = tensor::zero;
            }
            else
            {
                betaPatch[faceI] = 0;
                curDEpsEPred = DEpsilonPatch[faceI];
 
                if
                ( 
                    (DEpsilonEqPatch[faceI] >= 0)
                 && (sigmaEqEPatch[faceI] >= sigmaYPatch[faceI]) 
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
                               /(2*muPatch[faceI] + 3*lambdaPatch[faceI] + SMALL)
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
 
    // Update plastic strain increment
    scalar rf = 
        readScalar(plasticityModelCoeffs_.lookup("relaxationFactor"));

    volSymmTensorField newDEpsilonP = 
        4.5*beta_*mu_*(devSigma && DEpsilon)*devSigma
       /(
            (Ep_ + 3*mu_)*sqr(sigmaEq)
          + dimensioned<scalar>
            (
                "SMALL",
                mu_.dimensions()*sigmaEq.dimensions()*sigmaEq.dimensions(),
                SMALL
            )
        );

    DEpsilonP_ = rf*newDEpsilonP + (1.0 - rf)*DEpsilonP_;

    DEpsilonP_.correctBoundaryConditions();

    Info << "done" << endl;
}


void plasticityModel::updateYieldStress()
{
    Info << "Updating yield stress ... ";
/*
    const volScalarField epsilonEq =
        sqrt((2.0/3.0)*magSqr(dev(epsilon_ )))
      + dimensionedScalar("SMALL", dimless, SMALL);

    volScalarField Ep_ = Ep(epsilonEq);
*/
    const volSymmTensorField& newSigma = sigma();
    const volScalarField sigmaEq = sqrt(1.5*magSqr(dev(newSigma)));

    volScalarField Ep_ = Ep(sigmaEq);

    const scalarField& EpI = Ep_.internalField();
    const scalarField& sigmaEqI = sigmaEq.internalField();
    scalarField& sigmaYI = sigmaY_.internalField();

    forAll(sigmaYI, cellI)
    {
        if(EpI[cellI] != 0)
        {
            if( sigmaEqI[cellI] > sigmaYI[cellI] )
            {
                sigmaYI[cellI] = sigmaEqI[cellI];

                Info << " Internal cell " << cellI 
                    << " Yield stress updated to Sy= "
                    << sigmaEqI[cellI] * 1.0E-06 << " MPa" 
                    << endl;
            }
        }
    }

    forAll(sigmaY_.boundaryField(), patchI)
    {
        if (!sigmaY_.boundaryField()[patchI].coupled())
	{
        const scalarField& EpPatch = Ep_.boundaryField()[patchI];
        const scalarField& sigmaEqPatch = sigmaEq.boundaryField()[patchI];
        scalarField& sigmaYPatch = sigmaY_.boundaryField()[patchI];

        forAll(sigmaYPatch, faceI)
        {
            if(EpPatch[faceI] != 0)
            {
                if(sigmaEqPatch[faceI] > sigmaYPatch[faceI])
                {
                    sigmaYPatch[faceI] = sigmaEqPatch[faceI];

                    Info << "Boundary cell " << patchI << " " << faceI
                        << " Yield stress updated to Sy= "
                        << sigmaEqPatch[faceI] * 1.0E-06 << " MPa" 
                        << endl;
                }	
            }
        }
	}
    }	

    Info << "done" << endl;
}

bool plasticityModel::read()
{

    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
