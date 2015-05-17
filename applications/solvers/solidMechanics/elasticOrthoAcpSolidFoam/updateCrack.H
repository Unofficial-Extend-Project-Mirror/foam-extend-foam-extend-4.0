nFacesToBreak = 0;
nCoupledFacesToBreak = 0;
{
    // Check internal faces

    // scalarField effTraction =
    //   cohesiveZone.internalField() *
    //   mag(traction.internalField());
    scalarField normalTraction =
        cohesiveZone.internalField()*
        ( n.internalField() & traction.internalField() );

    // only consider tensile tractions
    normalTraction = max(normalTraction, scalar(0));
    scalarField shearTraction =
        cohesiveZone.internalField()*
        mag( (I - Foam::sqr(n.internalField())) & traction.internalField() );

    // the traction fraction is monitored to decide which faces to break:
    // ie (tN/tNC)^2 + (tS/tSC)^2 >1 to crack a face

    const surfaceScalarField sigmaMax = rheology.cohLaw().sigmaMax();
    const surfaceScalarField tauMax = rheology.cohLaw().tauMax();

    const scalarField& sigmaMaxI = sigmaMax.internalField();
    const scalarField& tauMaxI = tauMax.internalField();

    //scalarField effTractionFraction = effTraction/sigmaMax;
    scalarField effTractionFraction(normalTraction.size(), 0.0);

    if(cohesivePatchUPtr)
    {
        effTractionFraction =
            (normalTraction/sigmaMaxI)*(normalTraction/sigmaMaxI)
          + (shearTraction/tauMaxI)*(shearTraction/tauMaxI);
    }
    else
    {
        // solidCohesiveFixedModeMix only uses sigmaMax
        effTractionFraction =
            (normalTraction/sigmaMaxI)*(normalTraction/sigmaMaxI)
          + (shearTraction/sigmaMaxI)*(shearTraction/sigmaMaxI);
    }

    maxEffTractionFraction = gMax(effTractionFraction);

    SLList<label> facesToBreakList;
    SLList<scalar> facesToBreakEffTractionFractionList;

    forAll(effTractionFraction, faceI)
    {
        if (effTractionFraction[faceI] > 1.0)
        {
            facesToBreakList.insert(faceI);
            facesToBreakEffTractionFractionList.insert
            (
                effTractionFraction[faceI]
            );
        }
    }

    labelList facesToBreak(facesToBreakList);
    List<scalar> facesToBreakEffTractionFraction
    (
        facesToBreakEffTractionFractionList
    );

    nFacesToBreak = facesToBreak.size();

    // Break only one face per topo change
    if (nFacesToBreak > 1)
    {
        nFacesToBreak = 1;
    }

    // philipc - select face with maximum effective traction fraction
    label faceToBreakIndex = -1;
    scalar faceToBreakEffTractionFraction = 0;
    forAll(facesToBreakEffTractionFraction, faceI)
    {
        if
        (
            facesToBreakEffTractionFraction[faceI]
          > faceToBreakEffTractionFraction
        )
        {
            faceToBreakEffTractionFraction =
                facesToBreakEffTractionFraction[faceI];
            faceToBreakIndex = facesToBreak[faceI];
        }
    }

    scalar gMaxEffTractionFraction =
        returnReduce(faceToBreakEffTractionFraction, maxOp<scalar>());

    if (Pstream::parRun())
    {
        bool procHasFaceToBreak = false;
        if (nFacesToBreak > 0)
        {
            if
            (
                mag(gMaxEffTractionFraction - faceToBreakEffTractionFraction)
              < SMALL
            )
            {
                // philipc - Maximum traction fraction is on this processor
                procHasFaceToBreak = true;
            }
        }

        // Check if maximum is present on more then one processors
        label procID = Pstream::nProcs();
        if (procHasFaceToBreak)
        {
            procID = Pstream::myProcNo();
        }

        label minProcID =
            returnReduce<label>(procID, minOp<label>());

        if (procID != minProcID)
        {
            nFacesToBreak = 0;
        }
    }


    // Check coupled (processor) patches

    SLList<label> coupledFacesToBreakList;
    SLList<scalar> coupledFacesToBreakEffTractionFractionList;
    forAll(mesh.boundary(), patchI)
    {
        if (mesh.boundary()[patchI].coupled())
        {
            // scalarField pEffTraction =
            //   cohesiveZone.boundaryField()[patchI]*
            //   mag(traction.boundaryField()[patchI]);
            // scalarField pEffTractionFraction = pEffTraction/sigmaMax.boundaryField()[patchI];

            scalarField pNormalTraction =
                cohesiveZone.boundaryField()[patchI]*
                ( n.boundaryField()[patchI] & traction.boundaryField()[patchI] );

            // only consider tensile tractions
            pNormalTraction = max(pNormalTraction, scalar(0));

            scalarField pShearTraction =
                cohesiveZone.boundaryField()[patchI]*
                mag( (I - Foam::sqr(n.boundaryField()[patchI])) & traction.boundaryField()[patchI] );

            // the traction fraction is monitored to decide which faces to break:
            // ie (tN/tNC)^2 + (tS/tSC)^2 >1 to crack a face
            const scalarField& pSigmaMax = sigmaMax.boundaryField()[patchI];
            const scalarField& pTauMax = tauMax.boundaryField()[patchI];

            scalarField pEffTractionFraction(pNormalTraction.size(), 0.0);
            if(cohesivePatchUPtr)
            {
                pEffTractionFraction =
                    (pNormalTraction/pSigmaMax)*(pNormalTraction/pSigmaMax)
                  + (pShearTraction/pTauMax)*(pShearTraction/pTauMax);
            }
            else
            {
                // solidCohesiveFixedModeMix only uses sigmaMax
                pEffTractionFraction =
                    (pNormalTraction/pSigmaMax)*(pNormalTraction/pSigmaMax)
                  + (pShearTraction/pSigmaMax)*(pShearTraction/pSigmaMax);
            }

            label start = mesh.boundaryMesh()[patchI].start();

            forAll(pEffTractionFraction, faceI)
            {
                if (pEffTractionFraction[faceI] > maxEffTractionFraction)
                {
                    maxEffTractionFraction = pEffTractionFraction[faceI];
                }

                if (pEffTractionFraction[faceI] > 1.0)
                {
                    coupledFacesToBreakList.insert(start + faceI);
                    coupledFacesToBreakEffTractionFractionList.insert
                    (
                        pEffTractionFraction[faceI]
                    );
                }
            }
        }
    }

    labelList coupledFacesToBreak(coupledFacesToBreakList);
    List<scalar> coupledFacesToBreakEffTractionFraction
    (
        coupledFacesToBreakEffTractionFractionList
    );

    nCoupledFacesToBreak = coupledFacesToBreak.size();

    // Break only one face per topo change
    if (nCoupledFacesToBreak > 1)
    {
        nCoupledFacesToBreak = 1;
    }

    // Select coupled face with maximum effective traction fraction
    label coupledFaceToBreakIndex = -1;
    scalar coupledFaceToBreakEffTractionFraction = 0;
    forAll(coupledFacesToBreakEffTractionFraction, faceI)
    {
        if
        (
            coupledFacesToBreakEffTractionFraction[faceI]
          > coupledFaceToBreakEffTractionFraction
        )
        {
            coupledFaceToBreakEffTractionFraction =
                coupledFacesToBreakEffTractionFraction[faceI];
            coupledFaceToBreakIndex = coupledFacesToBreak[faceI];
        }
    }

    scalar gMaxCoupledEffTractionFraction =
        returnReduce(coupledFaceToBreakEffTractionFraction, maxOp<scalar>());

    if (Pstream::parRun())
    {

        bool procHasCoupledFaceToBreak = false;
        if (nCoupledFacesToBreak > 0)
        {
            if
            (
                mag(gMaxCoupledEffTractionFraction - coupledFaceToBreakEffTractionFraction)
              < SMALL
            )
            {
                // Maximum traction fraction is on this processor
                procHasCoupledFaceToBreak = true;
            }
        }

        // Check if maximum is present on more then one processors
        label procID = Pstream::nProcs();
        if (procHasCoupledFaceToBreak)
        {
            procID = Pstream::myProcNo();
        }

        label minProcID =
            returnReduce<label>(procID, minOp<label>());

        if (procID != minProcID)
        {
            nCoupledFacesToBreak = 0;
        }
    }

    if (gMaxCoupledEffTractionFraction > gMaxEffTractionFraction)
    {
        // Break coupled face
        nFacesToBreak = 0;
    }
    else
    {
        // Break internal face
        nCoupledFacesToBreak = 0;
    }

    // Make sure that coupled faces are broken in pairs
    labelList ngbProc(Pstream::nProcs(), -1);
    labelList index(Pstream::nProcs(), -1);
    if (nCoupledFacesToBreak)
    {
        label patchID =
            mesh.boundaryMesh().whichPatch(coupledFaceToBreakIndex);

        label start = mesh.boundaryMesh()[patchID].start();
        label localIndex = coupledFaceToBreakIndex - start;

        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(mesh.boundaryMesh()[patchID]);
        label ngbProcNo = procPatch.neighbProcNo();

        ngbProc[Pstream::myProcNo()] = ngbProcNo;
        index[Pstream::myProcNo()] = localIndex;
    }

    if (returnReduce(nCoupledFacesToBreak, maxOp<label>()))
    {
        reduce(ngbProc, maxOp<labelList>());
        reduce(index, maxOp<labelList>());

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                if (ngbProc[procI] == Pstream::myProcNo())
                {
                    forAll(mesh.boundaryMesh(), patchI)
                    {
                        if
                        (
                            mesh.boundaryMesh()[patchI].type()
                         == processorPolyPatch::typeName
                        )
                        {
                            const processorPolyPatch& procPatch =
                                refCast<const processorPolyPatch>
                                (
                                    mesh.boundaryMesh()[patchI]
                                );
                            label ngbProcNo = procPatch.neighbProcNo();

                            if (ngbProcNo == procI)
                            {
                                label start =
                                    mesh.boundaryMesh()[patchI].start();
                                coupledFaceToBreakIndex = start + index[procI];
                                nCoupledFacesToBreak = 1;
                            }
                        }
                    }
                }
            }
        }
    }


    vector faceToBreakTraction = vector::zero;
    vector faceToBreakNormal = vector::zero;
    scalar faceToBreakSigmaMax = 0.0;
    scalar faceToBreakTauMax = 0.0;

    // Set faces to break
    if (nFacesToBreak > 0)
    {
        faceToBreakTraction = traction.internalField()[faceToBreakIndex];
        faceToBreakNormal = n.internalField()[faceToBreakIndex];

        // Scale broken face traction
        faceToBreakSigmaMax = sigmaMaxI[faceToBreakIndex];
        faceToBreakTauMax = tauMaxI[faceToBreakIndex];
        scalar normalTrac = faceToBreakNormal & faceToBreakTraction;
        normalTrac = max(normalTrac, 0.0);
        scalar shearTrac = mag( (I - sqr(faceToBreakNormal)) & faceToBreakTraction );
        scalar scaleFactor = 1;
        if(cohesivePatchUPtr)
        {
            scaleFactor =
                Foam::sqrt
                (
                    1 /
                    (
                        (normalTrac/faceToBreakSigmaMax)*(normalTrac/faceToBreakSigmaMax)
                      + (shearTrac/faceToBreakTauMax)*(shearTrac/faceToBreakTauMax)
                    )
                );
        }
        else
        {
            // solidCohesiveFixedModeMix only uses sigmaMax
            scaleFactor =
                Foam::sqrt
                (
                    1 /
                    (
                        (normalTrac/faceToBreakSigmaMax)*(normalTrac/faceToBreakSigmaMax)
                      + (shearTrac/faceToBreakSigmaMax)*(shearTrac/faceToBreakSigmaMax)
                    )
                );
        }

        faceToBreakTraction *= scaleFactor;

        topoChange = true;
    }
    else if (nCoupledFacesToBreak > 0)
    {
        label patchID =
            mesh.boundaryMesh().whichPatch(coupledFaceToBreakIndex);
        label start = mesh.boundaryMesh()[patchID].start();
        label localIndex = coupledFaceToBreakIndex - start;

        faceToBreakTraction = traction.boundaryField()[patchID][localIndex];
        faceToBreakNormal = n.boundaryField()[patchID][localIndex];

        // Scale broken face traction
        faceToBreakSigmaMax = sigmaMax.boundaryField()[patchID][localIndex];
        faceToBreakTauMax = tauMax.boundaryField()[patchID][localIndex];
        scalar normalTrac = faceToBreakNormal & faceToBreakTraction;
        normalTrac = max(normalTrac, 0.0);
        scalar shearTrac = mag( (I - sqr(faceToBreakNormal)) & faceToBreakTraction );
        scalar scaleFactor = 1;
        if(cohesivePatchUPtr)
        {
            scaleFactor =
                Foam::sqrt
                (
                    1 /
                    (
                        (normalTrac/faceToBreakSigmaMax)*(normalTrac/faceToBreakSigmaMax)
                      + (shearTrac/faceToBreakTauMax)*(shearTrac/faceToBreakTauMax)
                    )
                );
        }
        else
        {
            // solidCohesiveFixedModeMix only uses sigmaMax
            scaleFactor =
                Foam::sqrt
                (
                    1 /
                    (
                        (normalTrac/faceToBreakSigmaMax)*(normalTrac/faceToBreakSigmaMax)
                      + (shearTrac/faceToBreakSigmaMax)*(shearTrac/faceToBreakSigmaMax)
                    )
                );
        }

        faceToBreakTraction *= scaleFactor;

        topoChange = true;
    }

    reduce(topoChange, orOp<bool>());

    labelList faceToBreak(nFacesToBreak, faceToBreakIndex);
    boolList faceToBreakFlip(nFacesToBreak, false);
    labelList coupledFaceToBreak
    (
        nCoupledFacesToBreak,
        coupledFaceToBreakIndex
    );

    reduce(nFacesToBreak, maxOp<label>());
    reduce(nCoupledFacesToBreak, maxOp<label>());

    if (nFacesToBreak || nCoupledFacesToBreak)
    {
        Pout << "Internal face to break: " << faceToBreak << endl;
        Pout << "Coupled face to break: " << coupledFaceToBreak << endl;

        mesh.setBreak(faceToBreak, faceToBreakFlip, coupledFaceToBreak);
        mesh.update();

        const labelList& faceMap = mesh.topoChangeMap().faceMap();
        label start = mesh.boundaryMesh()[cohesivePatchID].start();

        // mu = rheology.mu();
        // lambda = rheology.lambda();
        // muf = fvc::interpolate(mu);
        // lambdaf = fvc::interpolate(lambda);
        C = rheology.C();
        K = rheology.K();
        Cf = fvc::interpolate(C);
        Kf = fvc::interpolate(K);

        // we need to modify propertiess after cracking otherwise momentum equation is wrong
        // but solidInterface seems to hold some information about old mesh
        // so we will delete it and make another
        // we could probably add a public clearout function
        // create new solidInterface
        //Pout << "Creating new solidInterface" << endl;
        //delete solidInterfacePtr;
        //solidInterfacePtr = new solidInterface(mesh, rheology);
        // delete demand driven data as the mesh has changed
        if(rheology.solidInterfaceActive())
        {
            rheology.solInterface().clearOut();
            solidInterfacePtr->modifyProperties(Cf, Kf);
        }

        // Local crack displacement
        vectorField UpI =
            U.boundaryField()[cohesivePatchID].patchInternalField();
        vectorField oldUpI =
            U.oldTime().boundaryField()[cohesivePatchID].patchInternalField();

        // Global crack displacement
        vectorField globalUpI = mesh.globalCrackField(UpI);
        vectorField globalOldUpI = mesh.globalCrackField(oldUpI);

        // C and K field on new crack faces must be updated
        symmTensor4thOrderField CPI = C.boundaryField()[cohesivePatchID].patchInternalField();
        diagTensorField KPI = K.boundaryField()[cohesivePatchID].patchInternalField();
        symmTensor4thOrderField globalCPI = mesh.globalCrackField(CPI);
        diagTensorField globalKPI = mesh.globalCrackField(KPI);

        // cohesivePatchU.size()
        int cohesivePatchSize(cohesivePatchUPtr ? cohesivePatchUPtr->size() : cohesivePatchUFixedModePtr->size());

        // Initialise U for new cohesive face
        const labelList& gcfa = mesh.globalCrackFaceAddressing();
        label globalIndex = mesh.localCrackStart();
        // for (label i=0; i<cohesivePatchU.size(); i++)
        for (label i=0; i<cohesivePatchSize; i++)
        {
            label oldFaceIndex = faceMap[start+i];

            // If new face
            if (oldFaceIndex == faceToBreakIndex)
            {
                U.boundaryField()[cohesivePatchID][i] =
                    0.5
                   *(
                       globalUpI[globalIndex]
                     + globalUpI[gcfa[globalIndex]]
                    );
                U.oldTime().boundaryField()[cohesivePatchID][i] =
                    0.5
                   *(
                       globalOldUpI[globalIndex]
                     + globalOldUpI[gcfa[globalIndex]]
                    );

                // initialise C and K on new faces
                // set new face value to value of internal cell
                Cf.boundaryField()[cohesivePatchID][i] = globalCPI[globalIndex];
                Kf.boundaryField()[cohesivePatchID][i] = globalKPI[globalIndex];

                globalIndex++;
            }
            else
            {
                globalIndex++;
            }
        }

        // we must calculate grad using interface
        // U at the interface has not been calculated yet as interface.correct()
        // has not been called yet
        // not really a problem as gradU is correct in second outer iteration
        // as long as this does not cause convergence problems for the first iterations.
        // we should be able to calculate the interface displacements without
        // having to call interface.correct()
        // todo: add calculateInterfaceU() function
        // interface grad uses Gauss, we need least squares
        //gradU = solidInterfacePtr->grad(U);
        gradU = fvc::grad(U); // leastSquaresSolidInterface grad scheme
        //snGradU = fvc::snGrad(U);

#       include "calculateTraction.H"
        //if (nFacesToBreak || nCoupledFacesToBreak) mesh.write(); traction.write();

        // Initialise initiation traction for new cohesive patch face
        // for (label i=0; i<cohesivePatchU.size(); i++)
        for (label i=0; i<cohesivePatchSize; i++)
        {
            label oldFaceIndex = faceMap[start+i];

            // If new face
            if
            (
                (oldFaceIndex == faceToBreakIndex)
             || (oldFaceIndex == coupledFaceToBreakIndex)
            )
            {
                vector n0 =
                    mesh.Sf().boundaryField()[cohesivePatchID][i]
                    /mesh.magSf().boundaryField()[cohesivePatchID][i];
                //vector n1 = -n0;

                if ((n0 & faceToBreakNormal) > SMALL)
                {
                    traction.boundaryField()[cohesivePatchID][i] =
                        faceToBreakTraction;

                    traction.oldTime().boundaryField()[cohesivePatchID][i] =
                        faceToBreakTraction;

                    if(cohesivePatchUPtr)
                    {
                        cohesivePatchUPtr->traction()[i] = faceToBreakTraction;
                    }
                    else
                    {
                        cohesivePatchUFixedModePtr->traction()[i] =
                            faceToBreakTraction;

                        cohesivePatchUFixedModePtr->initiationTraction()[i] =
                            faceToBreakTraction;
                    }
                }
                else
                {
                    traction.boundaryField()[cohesivePatchID][i] =
                        -faceToBreakTraction;
                    traction.oldTime().boundaryField()[cohesivePatchID][i] =
                        -faceToBreakTraction;

                    //cohesivePatchU.traction()[i] = -faceToBreakTraction;
                    if(cohesivePatchUPtr)
                    {
                        cohesivePatchUPtr->traction()[i] = -faceToBreakTraction;
                    }
                    else
                    {
                        cohesivePatchUFixedModePtr->traction()[i] =
                           -faceToBreakTraction;
                        cohesivePatchUFixedModePtr->initiationTraction()[i] =
                            -faceToBreakTraction;
                    }
                }
            }
        }

        // hmmnn we only need a reference for very small groups of cells
        // turn off for now
        //#       include "updateReference.H"
    }
}
