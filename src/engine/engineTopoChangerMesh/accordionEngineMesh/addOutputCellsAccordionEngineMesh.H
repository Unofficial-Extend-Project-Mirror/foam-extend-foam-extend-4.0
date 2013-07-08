    {
        word cylinderOutputName(engTime().engineDict().lookup("cylinderSetName"));

        IOobject outputHeader
        (
            cylinderOutputName,
            time().constant()+"/polyMesh/sets/",
            *this,
            IOobject::MUST_READ
        );

        if(!outputHeader.headerOk())
        {
            WarningIn
            (
                "accordionEngineMesh::addZonesAndModifiers()"
            )   << "faceSet called " << cylinderOutputName
                << " does not exist. Continuing mesh motion without cylinder output cells" << endl;

        cz.append
        (
          new cellZone
          (
                "cylinderOutputCells",
                labelList(0),
                nCellZones,
                cellZones()
          )
      );

      nCellZones++;


        }
        else
        {

        cellSet outputCellSet(*this, cylinderOutputName);

        labelList outputCellSetList = outputCellSet.toc();

        cz.append
        (
          new cellZone
          (
                "cylinderOutputCells",
                outputCellSetList,
                nCellZones,
                cellZones()
          )
      );

      nCellZones++;
    }
}
