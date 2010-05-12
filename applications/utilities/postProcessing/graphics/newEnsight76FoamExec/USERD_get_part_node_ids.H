int USERD_get_part_node_ids
(
    int part_number,
    int *nodeid_array
)
{
#ifdef ENSIGHTDEBUG
    Info << "Entering: USERD_get_part_node_ids" << endl 
        << "part_number = " << part_number << endl 
        << flush;
#endif

    if (part_number == 1)
    {
        label nPoints = meshPtr->points().size();

        for (label indx=0; indx<nPoints; indx++)
        {
            nodeid_array[indx] = indx + 1;
        }
    }
    else if  (part_number < nPatches+2)
    {

        label patchi = part_number-2;
        const polyBoundaryMesh& bMesh = meshPtr->boundaryMesh();
        const vectorField& points = bMesh[patchi].points();

        label nPoints = points.size();
        
        for (label indx=0; indx<nPoints; indx++)
        {
            nodeid_array[indx] = indx + 1;
        }

    }
    else if  (sprayPtr && part_number == nPatches+2)
    {
        label indx = 0;
        for
        (
            Cloud<passiveParticle>::iterator elmnt = sprayPtr->begin();
            elmnt != sprayPtr->end();
            ++elmnt
        )
        {
             nodeid_array[indx] = indx + 1;
             indx++;
        }
    }
    else if
    (
        secondMeshPtr
     && part_number == secondMeshPartNum
    )
    {
        label nPoints = secondMeshPtr->points().size();

        for (label indx=0; indx<nPoints; indx++)
        {
            nodeid_array[indx] = indx + 1;
        }
    }
    else if
    (
        secondMeshPtr
     && part_number <= secondMeshPartNum + nSecondMeshPatches
    )
    {
        label patchi = part_number - secondMeshPartNum - 1;

        const polyBoundaryMesh& bMesh = secondMeshPtr->boundaryMesh();
        const vectorField& points = bMesh[patchi].points();
        label nPoints = points.size();

        for (label indx=0; indx<nPoints; indx++)
        {
            nodeid_array[indx] = indx + 1;
        }
    }
    else if
    (
        faMeshPtr
     && part_number == Numparts_available
    )
    {
        label nPoints = faMeshPtr->points().size();

        for (label indx=0; indx<nPoints; indx++)
        {
            nodeid_array[indx] = indx + 1;
        }
    }
    else
    {
        return Z_ERR;
    }

#ifdef ENSIGHTDEBUG
    Info << "Leaving: USERD_get_part_node_ids" << endl 
        << flush;
#endif

    return Z_OK;
}
