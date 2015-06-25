/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description
    Reads the specified surface and writes it in the fms format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "fileName.H"
#include "triSurf.H"
#include "triSurfaceChecks.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("input surface file");
    argList args(argc, argv);

    const fileName inFileName(args.args()[1]);

    triSurf surf(inFileName);

    label nFailed(0);

    boundBox bb;
    triSurfaceChecks::calculateBoundingBox(surf, bb);
    Info << "\nNumber of points " << surf.nPoints() << endl;
    Info << "Number of triangles " << surf.size() << endl;
    Info << "Number of patches " << surf.patches().size() << endl;
    Info << "Number of feature edges " << surf.nFeatureEdges() << endl;
    Info << "Bounding box " << bb << nl << nl << endl;

    const scalar distTol = SMALL * bb.mag();

    //- calculate manifolds
    const label nManifolds = triSurfaceChecks::checkSurfaceManifolds(surf);
    if( nManifolds > 1 )
    {
        ++nFailed;

        Info << "\nSurface mesh consists of " << nManifolds
             << " manifolds!!" << endl;
        Info << "You cannot mesh geometries consisting of more than"
                << " one domain, and it must not contain baffles"
                << " in the domain which shall be meshed." << endl;
    }
    else
    {
        Info << "\nSurface mesh consists of a single manifold." << endl;
    }

    //- find open boundary edges
    if( triSurfaceChecks::checkForHoles(surf) )
    {
        ++nFailed;

        Info << "\nSurface mesh has open boundaries!!" << endl;
        Info << "This indicates that there may be some holes in the surface"
             << " mesh. Holes in the mesh must be smaller than the specified"
             << " cell size at this location. In addition, please avoid"
             << " using the automatic refinement procedure."
             << " Please avoid using the minCellSize option." << endl;
    }
    else
    {
        Info << "No open edges found in the surface mesh." << endl;
    }

    //- find non-manifold edges
    if( triSurfaceChecks::checkForNonManifoldEdges(surf) )
    {
        ++nFailed;

        Info << "\nSurface mesh has non-manifold edges!!" << endl;
        Info << "This indicates that the surface mesh consists of multiple"
             << " domains and/or baffles. Please make sure that they are not"
             << " in the domain which shall be meshed." << endl;
    }
    else
    {
        Info << "Surface does not have any non-manifold edges." << endl;
    }

    //- check the number of disconnected parts
    if( triSurfaceChecks::checkDisconnectedParts(surf) > 1 )
    {
        ++nFailed;

        Info << "\nSurface mesh consists of disconnected parts!!" << endl;
        Info << "This is not a problem if there exists a region surrounding"
             << " all other regions! In other case, the mesher will generate"
             << " the mesh in the domains with most cells." << endl;
    }
    else
    {
        Info << "Surface mesh consists of a single region." << endl;
    }

    //- find triangles with small angles
    if( triSurfaceChecks::checkAngles(surf, "smallAngles", 1.0) )
    {
        ++nFailed;

        Info << "\nSurface mesh has some bad-quality triangles with"
             << " angles smaller than 1.0 deg!!" << endl;
        Info << "This may cause problems to the automatic refinement"
             << " procedure. Please avoid using the minCellSize option."
             << endl;
    }
    else
    {
        Info << "No sliver triangles found." << endl;
    }

    //- find self-intersections in the surface mesh
    if
    (
        triSurfaceChecks::checkSelfIntersections
        (
            surf,
            "selfIntersect",
            distTol
        )
    )
    {
        ++nFailed;

        Info << "\nFound self-intersecting parts in the surface mesh!!" << endl;
        Info << "This causes problems to the automatic refinement procedure"
                << " Please avoid using the minCellSize option."
                << " It can also cause problems to the boundary layer"
                << " generation procedure." << endl;
    }
    else
    {
        Info << "No self-intersections found." << endl;
    }

    //- find overlaps in the surface mesh
    if
    (
        triSurfaceChecks::checkOverlaps
        (
            surf,
            "overlaps",
            distTol,
            5.0
        )
    )
    {
        ++nFailed;

        Info << "\nFound overlapping parts in the surface mesh!!" << endl;
        Info << "This causes problems to the automatic refinement procedure."
             << " Please avoid using the minCellSize option." << endl;
    }

    //- check for existence of collocated points
    if
    (
        triSurfaceChecks::checkCollocatedPoints
        (
            surf,
            "collocatedPoints",
            distTol
        )
    )
    {
        ++nFailed;

        Info << "\nFound collocated points in the surface mesh!!" << endl;
        Info << "This causes problems to the automatic refinement procedure."
             << " Please avoid using the minCellSize option." << endl;
    }

    Info << nl << endl;

    if( nFailed )
    {
        Info << "\nFound " << nFailed
                << " checks indicating potential problems." << endl;
        Info << "However, it does not mean that you cannot generate"
                << " a valid mesh.\n" << endl;

        surf.writeSurface(inFileName);
    }
    else
    {
        Info << "\nSurface passes all checks.\n" << endl;
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
