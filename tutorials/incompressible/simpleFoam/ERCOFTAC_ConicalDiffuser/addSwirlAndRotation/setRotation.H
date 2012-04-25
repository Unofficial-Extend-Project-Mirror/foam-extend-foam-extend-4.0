// Add rotation/swirl at rotating patches
if (modifyBoundaries)
{
    forAll (rotPatches, patchi)
    {
        label patchID =
              mesh.boundaryMesh().findPatchID(rotPatches[patchi]);

        if (patchID > -1)
        {
            const fvPatch& cPatch = mesh.boundary()[patchID];
            const vectorField& faceCentres = cPatch.Cf();
            U.boundaryField()[patchID] ==
            (
	        U.boundaryField()[patchID]
                + (omega.value() ^ (faceCentres - center.value()))
            );
        }
        else
        {
        Info<< "Patch " << rotPatches[patchi] << " not found" << endl;
        }
    }
}

// Add swirl at internal cells
if (modifyInterior)
{
    forAll(mesh.C(), celli)
    {
        U[celli] += omega.value() ^ (mesh.C()[celli] - center.value());
    }
}
