    boolList needRef(mesh.regions().nRegions(), true);
    labelList refCells(mesh.regions().nRegions(), -1);

    labelList deadCells;

#   include "updateReference.H"

