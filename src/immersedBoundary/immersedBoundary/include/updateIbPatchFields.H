// Update all IB patch fields on mesh motion
// For moving mesh or moving IB surface, the intersection of IB with the
// background mesh changes, and IB patches are resized
// To allow smooth running of the code, IB patch fields need to be
// re-initialised on the new intersection
// HJ, 21/Jun/2018
{
    // Get reference to boundary
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Collect ib patch indices
    labelList ibPatches(patches.size());
    label nIbPatches = 0;

    forAll (patches, patchI)
    {
        if (isA<immersedBoundaryPolyPatch>(patches[patchI]))
        {
            ibPatches[nIbPatches] = patchI;
            nIbPatches++;
        }
    }

    // Resize ib patch index list
    ibPatches.setSize(nIbPatches);

    // Collect all vol and surface fields and evaluate the immersed boundary
    // patch field condition if a resize is needed

    // Exclude named volScalarFields
    HashSet<word> excludeVsf;

    // Exclude named volScalarFields
    HashSet<word> excludeVvf;

    // Exclude named surfaceScalarFields
    HashSet<word> excludeSsf;
    excludeSsf.insert("weightingFactors");
    excludeSsf.insert("deltaCoeffs");
    excludeSsf.insert("meshPhi");

    // volScalarField
    HashTable<volScalarField*> vsf = mesh.lookupClass<volScalarField>();

    for
    (
        HashTable<volScalarField*>::iterator fieldIter = vsf.begin();
        fieldIter != vsf.end();
        ++fieldIter
    )
    {
        volScalarField& f = *fieldIter();

        if (!excludeVsf.found(f.name()))
        {
            if (immersedBoundaryFvPatch::debug)
            {
                Info<< "Updating vs field " << f.name() << ": ";
            }

            forAll (ibPatches, ibpI)
            {
                if
                (
                    f.boundaryField()[ibPatches[ibpI]].size()
                 != mesh.boundary()[ibPatches[ibpI]].size()
                )
                {
                    if (immersedBoundaryFvPatch::debug)
                    {
                        Info<< ".";
                    }

                    f.boundaryField()[ibPatches[ibpI]].updateOnMotion();
                }
            }

            if (immersedBoundaryFvPatch::debug)
            {
                Info<< endl;
            }
        }
    }

    // volVectorField
    HashTable<volVectorField*> vvf = mesh.lookupClass<volVectorField>();

    for
    (
        HashTable<volVectorField*>::iterator fieldIter = vvf.begin();
        fieldIter != vvf.end();
        ++fieldIter
    )
    {
        volVectorField& f = *fieldIter();

        if (!excludeVvf.found(f.name()))
        {
            if (immersedBoundaryFvPatch::debug)
            {
                Info<< "Updating vv field " << f.name() << ": ";
            }

            forAll (ibPatches, ibpI)
            {
                if
                (
                    f.boundaryField()[ibPatches[ibpI]].size()
                 != mesh.boundary()[ibPatches[ibpI]].size()
                )
                {
                    if (immersedBoundaryFvPatch::debug)
                    {
                        Info<< ".";
                    }

                    f.boundaryField()[ibPatches[ibpI]].updateOnMotion();
                }
            }

            if (immersedBoundaryFvPatch::debug)
            {
                Info<< endl;
            }
        }
    }

    // surfaceScalarField
    HashTable<surfaceScalarField*> ssf = mesh.lookupClass<surfaceScalarField>();

    for
    (
        HashTable<surfaceScalarField*>::iterator fieldIter = ssf.begin();
        fieldIter != ssf.end();
        ++fieldIter
    )
    {
        surfaceScalarField& f = *fieldIter();

        if (!excludeSsf.found(f.name()))
        {
            if (immersedBoundaryFvPatch::debug)
            {
                Info<< "Updating ss field " << f.name() << ": ";
            }

            forAll (ibPatches, ibpI)
            {
                if
                (
                    f.boundaryField()[ibPatches[ibpI]].size()
                 != mesh.boundary()[ibPatches[ibpI]].size()
                )
                {
                    if (immersedBoundaryFvPatch::debug)
                    {
                        Info<< ".";
                    }

                    f.boundaryField()[ibPatches[ibpI]].updateOnMotion();
                }
            }

            if (immersedBoundaryFvPatch::debug)
            {
                Info<< endl;
            }
        }
    }
}
