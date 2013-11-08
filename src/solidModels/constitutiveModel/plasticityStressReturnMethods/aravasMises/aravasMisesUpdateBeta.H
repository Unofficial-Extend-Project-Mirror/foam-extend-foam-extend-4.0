    // Update internal beta
volScalarField& beta_ = activeYield_;

    scalarField& sigmaYI = sigmaY_.internalField();
    const volScalarField DEpsilonEq =
        sqrt((2.0/3.0)*magSqr(dev(newEpsilon)))
      + dimensionedScalar("SMALL", dimless, SMALL);
    // const volSymmTensorField& oldSigma = sigma();
    const volSymmTensorField& oldSigma =
      mesh.objectRegistry::lookupObject<volSymmTensorField>("sigma");
    const volSymmTensorField DSigma =
        2*mu*(DEpsilon - DEpsilonP_) + I*(lambda*tr(DEpsilon));
    // const volSymmTensorField sigma_ = sigma() + DSigma;
    const volSymmTensorField sigma = oldSigma + DSigma;
    volScalarField oldSigmaEq = sqrt(1.5*magSqr(dev(oldSigma)));
    const volSymmTensorField devSigma = dev(sigma);
    const volSymmTensorField DSigmaE = DSigma + 2*mu*DEpsilonP_;
    const volScalarField sigmaEqE  = sqrt(1.5*magSqr(dev(oldSigma + DSigmaE)));
    const volScalarField DSigmaEqE = sqrt(1.5*magSqr(dev(DSigmaE)));

    const scalarField& muI = mu.internalField();
    const scalarField& lambdaI = lambda.internalField();
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

        if ( (DEpsilonEqI[cellI] >= 0) && (oldBetaI[cellI] > SMALL) )
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

          // philipc -  use new for elasto-plastic steps
          // blending of new and old is only used for beta = 1 above
          //Ep_.boundaryField()[patchI][faceI] =
          //newEp.boundaryField()[patchI][faceI];
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
