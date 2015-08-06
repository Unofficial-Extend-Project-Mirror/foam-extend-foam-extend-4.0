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

Class
    triSurface2DCheck

Description


SourceFiles
    triSurface2DCheck.C

\*---------------------------------------------------------------------------*/

#ifndef triSurface2DCheck_H
#define triSurface2DCheck_H

#include "labelLongList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class triSurf;
class boundBox;

/*---------------------------------------------------------------------------*\
                    Namespace triSurfaceChecka Declaration
\*---------------------------------------------------------------------------*/

namespace triSurfaceChecks
{

//- check the angles of triangles. Triangles wth the smallest angle below
//- the given threshold are marked as bad-quality triangles
label checkAngles
(
    const triSurf&,
    labelLongList&,
    const scalar angleTol = 5.0
);

label checkAngles
(
    triSurf&,
    const word subsetName = "smallAngles",
    const scalar angleTol = 5.0
);

//- check for existence of boudary and non-manifold edges
label checkSurfaceManifolds(const triSurf&, labelLongList&triangleInManifold);
label checkSurfaceManifolds(triSurf&, const word subsetPrefix="manifold_");

//- check for existence of holes in the surface mesh
label checkForHoles(const triSurf&, labelLongList&);
label checkForHoles(triSurf&, const word subsetName="holes");

//- check for existence of non-manifold edges
label checkForNonManifoldEdges(const triSurf&, labelLongList&);
label checkForNonManifoldEdges
(
    triSurf&,
    const word subsetPrefix="nonManifold_"
);

//- check orientation of triangles
label checkOrientation(const triSurf&, labelLongList& triangleInGroup);
label checkOrientation(triSurf&, const word subsetPrefix="orientation_");

//- check for disconnected surface parts
label checkDisconnectedParts(const triSurf&, labelLongList& triangleInRegion);
label checkDisconnectedParts(triSurf&, const word subsetPrefix="region_");

//- calculate bounding box of the surface mesh
void calculateBoundingBox(const triSurf&, boundBox&);

//- check existence of collocated points
label checkCollocatedPoints
(
    const triSurf&,
    labelLongList&,
    const scalar distTol = 1e-6
);
label checkCollocatedPoints
(
    triSurf&,
    const word subsetName="collocatedPoints",
    const scalar distTol = 1e-6
);

//- check if there exist any self-intersections
label checkSelfIntersections
(
    const triSurf&,
    labelLongList&,
    const scalar tol=1e-6
);
label checkSelfIntersections
(
    triSurf&,
    const word subsetName="selfIntersect",
    const scalar tol=1e-6
);

//- check if there exist any overlaps
label checkOverlaps
(
    const triSurf&,
    labelLongList&,
    const scalar tol = 1e-6,
    const scalar angleTol = 5.0
);
label checkOverlaps
(
    triSurf&,
    const word subsetName="overlappingTriangles",
    const scalar tol = 1e-6,
    const scalar angleTol = 5.0
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace triSurfaceChecks

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
