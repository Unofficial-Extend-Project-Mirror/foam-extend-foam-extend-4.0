/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Create intermediate mesh files from SAMM files

\*---------------------------------------------------------------------------*/

#include "sammMesh.H"
#include "foamTime.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void sammMesh::writeMesh()
{
    if (isShapeMesh_)
    {
        Info << "This is a shapeMesh." << endl;

        polyMesh pShapeMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                runTime_.constant(),
                runTime_
            ),
            xferCopy(points_),           // we could probably re-use the data
            cellShapes_,
            boundary_,
            patchNames_,
            patchTypes_,
            defaultFacesName_,
            defaultFacesType_,
            patchPhysicalTypes_
        );

        Info << "Writing polyMesh" << endl;
        pShapeMesh.write();
    }
    else
    {
        // This is a polyMesh.

        createPolyMeshData();

        Info << "This is a polyMesh" << endl;

        polyMesh pMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                runTime_.constant(),
                runTime_
            ),
            xferCopy(points_),           // we could probably re-use the data
            xferCopy(meshFaces_),
            xferCopy(cellPolys_)
        );

        pMesh.addPatches(polyBoundaryPatches(pMesh));

        Info << "Writing polyMesh" << endl;
        pMesh.write();
    }
}


// ************************************************************************* //
