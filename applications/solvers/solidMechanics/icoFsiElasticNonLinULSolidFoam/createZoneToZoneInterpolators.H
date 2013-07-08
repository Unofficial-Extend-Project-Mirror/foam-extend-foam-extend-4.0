
zoneToZoneInterpolation* interpolatorFluidSolidPtr = NULL;
zoneToZoneInterpolation* interpolatorSolidFluidPtr = NULL;


if(!interpolatorFluidSolidPtr || !interpolatorSolidFluidPtr)
{
    deleteDemandDrivenData(interpolatorFluidSolidPtr);
    deleteDemandDrivenData(interpolatorSolidFluidPtr);

    Info << "Create fluid-to-solid and solid-to-fluid interpolators" << endl;

    interpolatorFluidSolidPtr = new zoneToZoneInterpolation
    (
        mesh.faceZones()[fluidZoneID](),
        stressMesh.faceZones()[solidZoneID](),
        intersection::VISIBLE
    );

    interpolatorSolidFluidPtr = new zoneToZoneInterpolation
    (
        stressMesh.faceZones()[solidZoneID](),
        mesh.faceZones()[fluidZoneID](),
        intersection::VISIBLE
    );

    Info << "Check fluid-to-solid and solid-to-fluid interpolators" << endl;

    {
        vectorField fluidPatchFaceCentres = 
            vectorField(mesh.boundaryMesh()[fluidPatchID].faceCentres());

        vectorField fluidZoneFaceCentres
        (
            mesh.faceZones()[fluidZoneID].size(), 
            vector::zero
        );

        
        const label fluidPatchStart = 
            mesh.boundaryMesh()[fluidPatchID].start();

        forAll (fluidPatchFaceCentres, i)
        {
            fluidZoneFaceCentres
            [
                mesh.faceZones()[fluidZoneID].whichFace(fluidPatchStart + i)
            ] = 
                fluidPatchFaceCentres[i];
        }

        // Parallel data exchange: collect faceCentres field on all processors
        reduce(fluidZoneFaceCentres, sumOp<vectorField>());

        vectorField solidZoneFaceCentres =
            interpolatorFluidSolidPtr->faceInterpolate
            (
                fluidZoneFaceCentres
            );

        vectorField solidPatchFaceCentres
            (
                stressMesh.boundaryMesh()[solidPatchID].size(), 
                vector::zero
            );

        const label solidPatchStart = 
            stressMesh.boundaryMesh()[solidPatchID].start();

        forAll(solidPatchFaceCentres, i)
        {
            solidPatchFaceCentres[i] =
                solidZoneFaceCentres
                [
                    stressMesh.faceZones()[solidZoneID]
                   .whichFace(solidPatchStart + i)
                ];
        }

        scalar maxDist = gMax
        (
            mag
            (
                solidPatchFaceCentres
              - stressMesh.boundaryMesh()[solidPatchID].faceCentres()
            )
        );

        Info << "Fluid-to-solid face interpolation error: " << maxDist 
            << endl;
    }



    {
        vectorField solidPatchFaceCentres = 
            vectorField(stressMesh.boundaryMesh()[solidPatchID].faceCentres());

        vectorField solidZoneFaceCentres
        (
            stressMesh.faceZones()[solidZoneID].size(), 
            vector::zero
        );

        const label solidPatchStart = 
            stressMesh.boundaryMesh()[solidPatchID].start();

        forAll (solidPatchFaceCentres, i)
        {
            solidZoneFaceCentres
            [
                stressMesh.faceZones()[solidZoneID]
               .whichFace(solidPatchStart + i)
            ] = 
                solidPatchFaceCentres[i];
        }

        // Parallel data exchange: collect faceCentres field on all processors
        reduce(solidZoneFaceCentres, sumOp<vectorField>());

        vectorField fluidZoneFaceCentres =
            interpolatorSolidFluidPtr->faceInterpolate
            (
                solidZoneFaceCentres
            );

        vectorField fluidPatchFaceCentres
        (
            mesh.boundaryMesh()[fluidPatchID].size(), 
            vector::zero
        );

        const label fluidPatchStart = 
            mesh.boundaryMesh()[fluidPatchID].start();

        forAll(fluidPatchFaceCentres, i)
        {
            fluidPatchFaceCentres[i] =
                fluidZoneFaceCentres
                [
                    mesh.faceZones()[fluidZoneID]
                   .whichFace(fluidPatchStart + i)
                ];
        }

        scalar maxDist = gMax
        (
            mag
            (
                fluidPatchFaceCentres
              - mesh.boundaryMesh()[fluidPatchID].faceCentres()
            )
        );

        Info << "Solid-to-fluid face interpolation error: " << maxDist
            << endl;
    }
}

zoneToZoneInterpolation& interpolatorFluidSolid = 
    *interpolatorFluidSolidPtr;

zoneToZoneInterpolation& interpolatorSolidFluid = 
    *interpolatorSolidFluidPtr;

