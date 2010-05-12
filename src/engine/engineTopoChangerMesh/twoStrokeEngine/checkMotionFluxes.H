{
    const fvsPatchScalarField& pistonMeshPhi =
        phi().boundaryField()[boundaryMesh().findPatchID("piston")];

    Info << "max piston mesh phi 2 = " << max(mag(pistonMeshPhi)) << endl;

    const labelList& cutFaces =
        faceZones()[faceZones().findZoneID("cutFaceZone")];

    forAll(cutFaces, facei)
    {
        correctedMeshPhi[cutFaces[facei]] = 0.0;
    }

//     setPhi(correctedMeshPhi);

    Info<< "max mesh phi, corrected cut faces = "
        << max(mag(phi().internalField())) << nl
        << "max mesh phi, corrected cut faces = "
        << min(mag(phi().internalField())) << endl;

    const polyPatch& wallPatch =
        boundaryMesh()[boundaryMesh().findPatchID("wall")];

    label wallSize = wallPatch.size();

    fvsPatchScalarField& wallMeshPhi =
        correctedMeshPhi.boundaryField()
        [
            boundaryMesh().findPatchID("wall")
        ];

    forAll(wallMeshPhi, k)
    {
        if(mag(wallMeshPhi[k]) > SMALL)
        {
            wallMeshPhi[k] = 0;
        }
    }

    fvsPatchScalarField& linerMeshPhi =
        correctedMeshPhi.boundaryField()[boundaryMesh().findPatchID("liner")];

    forAll(linerMeshPhi, k)
    {
        if(mag(linerMeshPhi[k]) > SMALL)
        {
            linerMeshPhi[k] = 0;
        }
    }

//        setPhi(correctedMeshPhi);

    Info<< "max mesh phi after correction = " << max(mag(phi())) << nl
        << "min mesh phi after correction = " << min(mag(phi())) << endl;

//    if(1 == 0)
    {
        scalarField dVdt = (1.0 - V0()/V())/engTime().deltaT();

        scalarField sumMeshPhi(V().size(), 0.0);
        scalarField sumMeshPhiCorrected(V().size(), 0.0);

        forAll(owner(), facei)
        {
            sumMeshPhi[owner()[facei]] += phi()[facei];
            sumMeshPhi[neighbour()[facei]] -= phi()[facei];
            sumMeshPhiCorrected[owner()[facei]] += correctedMeshPhi[facei];
            sumMeshPhiCorrected[neighbour()[facei]] -= correctedMeshPhi[facei];
        }

        forAll(boundary(), patchi)
        {
            const unallocLabelList& pFaceCells =
                boundary()[patchi].faceCells();

            const fvsPatchField<scalar>& pssf = phi().boundaryField()[patchi];
            const fvsPatchField<scalar>& pssfCorrected =
                correctedMeshPhi.boundaryField()[patchi];

            forAll(boundary()[patchi], facei)
            {
                sumMeshPhi[pFaceCells[facei]] += pssf[facei];
                sumMeshPhiCorrected[pFaceCells[facei]] += pssfCorrected[facei];
            }
        }

        sumMeshPhi /= V();

        sumMeshPhiCorrected /= V();

        scalarField sumCheck = dVdt - sumMeshPhi;
        scalarField sumCheckCorrected = dVdt - sumMeshPhiCorrected;

        label nOutCheck = 0;
        label nOutCheckCorrected = 0;

        forAll(sumCheck, celli)
        {
            if(mag(sumCheck[celli]) > 1)
            {
                nOutCheck++;
            }
        }

        forAll(sumCheckCorrected, celli)
        {
            if (mag(sumCheckCorrected[celli]) > 1)
            {
                nOutCheckCorrected++;
            }
        }

        Info<< "found " << nOutCheck
            << " cells with inconsistent motion fluxes real" << endl;
        Info<< "found " << nOutCheckCorrected
            << " cells with inconsistent motion fluxes CORRECTED" << endl;
    }

    /*
        newPoints = allPoints();

//        movePoints(newPoints);

//        pointField mappedOldPointsNew(newPoints.size());
        pointField mappedOldPointsNew(points().size());

        mappedOldPointsNew.map(oldPointsNew, topoChangeMap4->pointMap());

        // Solve the correct mesh motion to make sure motion fluxes
        // are solved for and not mapped
        movePoints(mappedOldPointsNew);

//        resetMotion();

        Info << "max mesh phi after reset motion = " << max(mag(phi())) << endl;
        Info << "max mesh phi after reset motion = " << min(mag(phi())) << endl;

        setV0();
        movePoints(newPoints);
//        pointField newPointsNew = allPoints();
//        movePoints(newPointsNew);

        Info << "max mesh phi after corr = " << max(mag(phi())) << endl;
        Info << "max mesh phi after corr = " << min(mag(phi())) << endl;

        */

        // Changing topology by hand
//        autoPtr<mapPolyMesh> topoChangeMap3 = topoChanger_.changeMesh();

    Info << "Sliding interfaces coupled: " << attached() << endl;
    Info << "max mesh phi after correction = " << max(mag(phi())) << endl;
}
