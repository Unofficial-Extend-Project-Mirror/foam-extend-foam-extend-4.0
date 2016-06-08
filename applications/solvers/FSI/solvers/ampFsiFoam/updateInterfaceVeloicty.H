{
    Info << max(mag(fsi.flow().U().boundaryField()[fsi.fluidPatchIndex()]))
        << ", "
        << average(mag(fsi.flow().U().boundaryField()[fsi.fluidPatchIndex()]))
        << endl;

    const volVectorField& cU = fsi.flow().U();
    volVectorField& U = const_cast<volVectorField& >(cU);
    // Store old time step velocity
    U.oldTime();

//     const fvPatchField<vector>& ciVelocity =
//         fsi.flow().U().boundaryField()[fsi.fluidPatchIndex()];
//     fvPatchField<vector>& iVelocity =
//         const_cast<fvPatchField<vector>& >(ciVelocity);

    fvPatchField<vector>& iVelocity =
           U.boundaryField()[fsi.fluidPatchIndex()];

    // Calculate velocity using position of face-centres
    // from previous and current postion and assuming Euler implicity scheme
    // Nomal component of velocity is calculated using mesh motion flux
    {
        const fvPatch& p = mesh.boundary()[fsi.fluidPatchIndex()];
        const polyPatch& pp = p.patch();
        const pointField& oldPoints = mesh.oldPoints();

        vectorField oldFc(pp.size());

        forAll(oldFc, i)
        {
            oldFc[i] = pp[i].centre(oldPoints);
        }

        // Get wall-parallel mesh motion velocity from geometry
        vectorField Up =
            (pp.faceCentres() - oldFc)/mesh.time().deltaT().value();

        // Mesh motion flux
        scalarField phip =
            p.patchField<surfaceScalarField, scalar>
            (
                fvc::meshPhi(fsi.flow().U())
            );

//         Up = iVelocity;

        vectorField n = p.nf();
        const scalarField& magSf = p.magSf();
        scalarField Un = phip/(magSf + VSMALL);

        // Adjust for surface-normal mesh motion flux
        iVelocity = (Up + n*(Un - (n & Up)));

//         iVelocity = Up;
    }

    Info << max(mag(fsi.flow().U().boundaryField()[fsi.fluidPatchIndex()]))
        << ", "
        << average(mag(fsi.flow().U().boundaryField()[fsi.fluidPatchIndex()]))
        << endl;
}
