    triSurface patch(prefix + ".stl");

    // Read patches
    List<Tuple2<word, label> > patchNames(patch.patches().size());

    forAll (patch.patches(), patchI)
    {
        patchNames[patchI] =
            Tuple2<word, label>(patch.patches()[patchI].name(), patchI);
    }

    Info << "Patches: " << patchNames << endl;

    // Create polyMesh
    faceList faces(patch.size());
    labelList faceLabels(patch.size());
    labelList faceRegion(patch.size());

    forAll (patch, faceI)
    {
        faces[faceI] = face(patch[faceI]);
        faceLabels[faceI] = faceI;
        faceRegion[faceI] = patch[faceI].region();
    }

    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        xferCopy(patch.points()),
        xferCopy(faces),
        xferCopy(labelList(0)),
        xferCopy(labelList(0)),
        false
    );

{
    // Add zones
    List<pointZone*> pz(0);
    List<faceZone*> fz(1);
    List<cellZone*> cz(0);

    fz[0] = new faceZone
    (
        "roof",
        faceLabels,
        boolList(patch.size(), false),
        0,
        mesh.faceZones()
    );

    mesh.addPatches(List<polyPatch*>(0), true);
    mesh.addZones(pz, fz, cz);
}
