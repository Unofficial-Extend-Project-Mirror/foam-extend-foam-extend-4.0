{
    // Multiple IB patches currently not supported.
    // This required IB geometry update to re-combine the corrected fields in
    // Geometry correction functions in immersedBoundaryFvPatch
    // HJ, 30/Nov/2017

    // Correct boundary patches to resize
    gamma.correctBoundaryConditions();
    sGamma.correctBoundaryConditions();

    scalarField& gammaIn = gamma.internalField();

    gammaIn = mesh.V()/mesh.cellVolumes();

    scalarField& sGammaIn = sGamma.internalField();

    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField magFaceAreas = mag(mesh.faceAreas());

    sGammaIn =  magSf.internalField()/
        scalarField::subField(magFaceAreas, mesh.nInternalFaces());

    forAll (mesh.boundary(), patchI)
    {
        if (!isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
        {
            sGamma.boundaryField()[patchI] =
                magSf.boundaryField()[patchI]/
                mesh.boundary()[patchI].patchSlice(magFaceAreas);

            gamma.boundaryField()[patchI] =
                sGamma.boundaryField()[patchI];
        }
        else
        {
            sGamma.boundaryField()[patchI] =
                scalarField(mesh.boundary()[patchI].size(), 1);

            gamma.boundaryField()[patchI] =
                scalarField(mesh.boundary()[patchI].size(), 1);
        }
    }

    // Adjust the mask in dead faces and dead cells
    forAll (mesh.boundary(), patchI)
    {
        if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
        {
            const immersedBoundaryFvPatch& ibPatch =
                refCast<const immersedBoundaryFvPatch>
                (
                    mesh.boundary()[patchI]
                );

            const labelList& deadCells = ibPatch.ibPolyPatch().deadCells();

            forAll (deadCells, dcI)
            {
                gammaIn[deadCells[dcI]] = 0;
            }

            const labelList& deadFaces = ibPatch.ibPolyPatch().deadFaces();

            forAll (deadFaces, dfI)
            {
                if (mesh.isInternalFace(deadFaces[dfI]))
                {
                    sGammaIn[deadFaces[dfI]] = 0;
                }
                else
                {
                    const label pI =
                        mesh.boundaryMesh().whichPatch(deadFaces[dfI]);

                    if (!isA<emptyFvPatch>(mesh.boundary()[pI]))
                    {
                        const label fI =
                            mesh.boundaryMesh()[pI].whichFace(deadFaces[dfI]);

                        sGamma.boundaryField()[pI][fI] = 0;
                    }
                }
            }
        }
    }
}
