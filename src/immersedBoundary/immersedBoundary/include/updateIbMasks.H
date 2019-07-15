{
    // Multiple IB patches currently not supported.
    // This required IB geometry update to re-combine the corrected fields in
    // Geometry correction functions in immersedBoundaryFvPatch
    // HJ, 30/Nov/2017

    // Correct boundary patches to resize
    gamma.correctBoundaryConditions();
    sGamma.correctBoundaryConditions();


    // Reset gamma and sGamma
    gamma == dimensionedScalar("one", dimless, 1);
    sGamma == dimensionedScalar("one", dimless, 1);

    scalarField& gammaIn = gamma.internalField();
    scalarField& sGammaIn = sGamma.internalField();

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
