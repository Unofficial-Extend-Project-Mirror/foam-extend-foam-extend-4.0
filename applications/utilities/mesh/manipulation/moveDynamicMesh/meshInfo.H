{
    Info<< "Mesh data: " << nl
        << "    total points  : " << mesh.allPoints().size() << nl
        << "    live points   : " << mesh.nPoints() << nl
        << "    total faces   : " << mesh.allFaces().size() << nl
        << "    live faces    : " << mesh.nFaces() << nl
        << "    internal faces: " << mesh.neighbour().size() << nl
        << "    cells         : " << mesh.nCells() << nl
        << "    patches       : " << mesh.boundaryMesh().size() << nl
        << "    point zones   : " << mesh.pointZones().size() << nl
        << "    face zones    : " << mesh.faceZones().size() << nl
        << "    cell zones    : " << mesh.cellZones().size() << endl;
}
