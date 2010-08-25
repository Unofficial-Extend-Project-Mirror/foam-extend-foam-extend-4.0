// Note: coord_array is 1-based.

int USERD_get_part_coords
(
    int part_number,
    float **coord_array
)
{

#ifdef ENSIGHTDEBUG
    Info << "Entering: USERD_get_part_coords" << endl << 
        "part_number = " << part_number << endl << flush;
#endif

    if (part_number == 1)
    {
        //#       include "checkForNewMesh.H"

        const vectorField& points = meshPtr->points();
        label nPoints = points.size();
        
        for (label indx=0; indx<nPoints; indx++)
        {
            coord_array[0][indx+1] = scalingFactor*(float)points[indx].x();
            coord_array[1][indx+1] = scalingFactor*(float)points[indx].y();
            coord_array[2][indx+1] = scalingFactor*(float)points[indx].z();
        }
    }
    else if (part_number < nPatches+2)
    {
        //#       include "checkForNewMesh.H"

        label patchi = part_number-2;
        const polyBoundaryMesh& bMesh = meshPtr->boundaryMesh();
        const vectorField& points = bMesh[patchi].points();
        label nPoints = points.size();
        
        for (label indx=0; indx<nPoints; indx++)
        {
            coord_array[0][indx+1] = scalingFactor*(float)points[indx].x();
            coord_array[1][indx+1] = scalingFactor*(float)points[indx].y();
            coord_array[2][indx+1] = scalingFactor*(float)points[indx].z();
        }

    }
    else if (sprayPtr && part_number == nPatches+2)
    {
        label indx = 1;

        for
        (
            Cloud<passiveParticle>::iterator elmnt = sprayPtr->begin();
            elmnt != sprayPtr->end();
            ++elmnt
        )
        {
            coord_array[0][indx] = scalingFactor*(float)elmnt().position().x();
            coord_array[1][indx] = scalingFactor*(float)elmnt().position().y();
            coord_array[2][indx] = scalingFactor*(float)elmnt().position().z();
            indx++;
        }
    }
    else if
    (
        secondMeshPtr
     && part_number == secondMeshPartNum 
    )
    {
        const vectorField& points = secondMeshPtr->points();
        label nPoints = points.size();
        
        for (label indx=0; indx<nPoints; indx++)
        {
            coord_array[0][indx+1] = scalingFactor*(float)points[indx].x();
            coord_array[1][indx+1] = scalingFactor*(float)points[indx].y();
            coord_array[2][indx+1] = scalingFactor*(float)points[indx].z();
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
            coord_array[0][indx+1] = scalingFactor*(float)points[indx].x();
            coord_array[1][indx+1] = scalingFactor*(float)points[indx].y();
            coord_array[2][indx+1] = scalingFactor*(float)points[indx].z();
        }
    }
    else if
    (
        faMeshPtr
     && part_number == Numparts_available
    )
    {
        const vectorField& points = faMeshPtr->points();
        label nPoints = points.size();
        
        for (label indx=0; indx<nPoints; indx++)
        {
            coord_array[0][indx+1] = scalingFactor*(float)points[indx].x();
            coord_array[1][indx+1] = scalingFactor*(float)points[indx].y();
            coord_array[2][indx+1] = scalingFactor*(float)points[indx].z();
        }
    }
    else
    {
        return Z_ERR;
    }

#ifdef ENSIGHTDEBUG
    Info << "Leaving: USERD_get_part_coords" << endl << flush;
#endif

    return Z_OK;
}
