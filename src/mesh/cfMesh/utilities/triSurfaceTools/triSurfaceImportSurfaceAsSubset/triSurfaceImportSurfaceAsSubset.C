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

\*---------------------------------------------------------------------------*/

#include "triSurfaceImportSurfaceAsSubset.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceImportSurfaceAsSubset::createOctree
(
    const triSurf& surf,
    meshOctree& octree
)
{

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfaceImportSurfaceAsSubset::triSurfaceImportSurfaceAsSubset(triSurf& surface)
:
    surf_(surface),
    octreePtr_(NULL)
{}

triSurfaceImportSurfaceAsSubset::~triSurfaceImportSurfaceAsSubset()
{
    deleteDemandDrivenData(octreePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceImportSurfaceAsSubset::addSurfaceAsSubset
(
    const triSurf& importSurf,
    const word& subsetName,
    const scalar angleTol
)
{
    if( !octreePtr_ )
    {
        octreePtr_ = new meshOctree(surf_);
        meshOctreeCreator(*octreePtr_).createOctreeWithRefinedBoundary
        (
            direction(20),
            15
        );
    }

    const pointField& points = surf_.points();
    const vectorField& fNornals = surf_.facetNormals();
    const vectorField& fCentres = surf_.facetCentres();

    labelList nearestTriangle(importSurf.size(), -1);

    //- check which triangles in the surface fit best to the centres of the
    //- triangles in the import surface
    const pointField& importSurfPoints = importSurf.points();
    const vectorField& importFaceCentres = importSurf.facetCentres();
    const vectorField& importFaceNormals = importSurf.facetNormals();
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(nearestTriangle, triI)
    {
        point np;
        scalar dSq;
        label nt, patch;

        octreePtr_->findNearestSurfacePoint
        (
            np,
            dSq,
            nt,
            patch,
            importFaceCentres[triI]
        );

        //- find the longest edge distance
        scalar maxEdgeDSq(0.);
        const labelledTri& tri = importSurf[triI];
        forAll(tri, pI)
        {
            const point& s = importSurfPoints[tri[pI]];
            const point& e = importSurfPoints[tri[(pI+1)%3]];

            maxEdgeDSq = max(maxEdgeDSq, magSqr(e - s));
        }

        //- check if the triangle has been found
        if( (nt < 0) || (dSq > 0.09 * maxEdgeDSq) )
        {
            Warning << "Could not find a matching triangle " << endl;
            Warning << "It seems that your surface meshes do not overlap" << endl;
            continue;
        }

        vector nTri = importFaceNormals[triI];
        const scalar magSqrTri = magSqr(nTri);

        //- skip sliver triangles
        if( magSqrTri < VSMALL )
            continue;

        vector normal = fNornals[nt];
        const scalar dSqNormal = magSqr(normal);

        //- skip sliver triangles
        if( dSqNormal < VSMALL )
            continue;

        if( ((nTri & normal) / (magSqrTri * dSqNormal)) > angleTol )
            nearestTriangle[triI] = nt;
    }

    meshOctree otherSurfOctree(importSurf);
    meshOctreeCreator(otherSurfOctree).createOctreeWithRefinedBoundary(20, 15);

    //- search for nearest facets in the import surface
    DynList<label> containedTriangles;
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40) private(containedTriangles)
    # endif
    forAll(surf_, triI)
    {
        //- find the bounding box and the ize of the triangle
        boundBox bb(fCentres[triI], fCentres[triI]);

        scalar maxEdgeDSq(0.);
        const labelledTri& tri = surf_[triI];
        forAll(tri, pI)
        {
            //- bounding box of the surface triangle
            bb.min() = min(bb.min(), points[tri[pI]]);
            bb.max() = max(bb.max(), points[tri[pI]]);

            const point& s = points[tri[pI]];
            const point& e = points[tri[(pI+1)%3]];

            maxEdgeDSq = max(maxEdgeDSq, magSqr(e - s));
        }

        //- find the nearest triangle in the surface which shall be imported
        otherSurfOctree.findTrianglesInBox(bb, containedTriangles);

        label nt(-1);
        scalar dSq(VGREAT);
        forAll(containedTriangles, ctI)
        {
            const point p =
                help::nearestPointOnTheTriangle
                (
                    containedTriangles[ctI],
                    importSurf,
                    fCentres[triI]
                );

            const scalar distSq  = magSqr(p - fCentres[triI]);

            if( distSq < dSq )
            {
                nt = containedTriangles[ctI];
                dSq = distSq;
            }
        }

        //- check if the triangle has been found
        if( (nt < 0) || (dSq > 0.09 * maxEdgeDSq) )
            continue;

        //- skip firther checkes f it has found the same triangle
        if( nearestTriangle[nt] == triI )
            continue;

        vector nTri = fNornals[triI];
        const scalar magSqrTri = magSqr(nTri);

        //- skip sliver triangles
        if( magSqrTri < VSMALL )
            continue;

        vector normal = importFaceNormals[nt];
        const scalar dSqNormal = magSqr(normal);

        //- skip sliver triangles
        if( dSqNormal < VSMALL )
            continue;

        if( ((nTri & normal) / (magSqrTri * dSqNormal)) > angleTol )
            nearestTriangle[nt] = triI;
    }

    //- create a facet subset in the surface mesh and add the facets into it
    const label subsetId = surf_.addFacetSubset(subsetName);

    forAll(nearestTriangle, triI)
    {
        if( nearestTriangle[triI] < 0 )
            continue;

        surf_.addFacetToSubset(subsetId, nearestTriangle[triI]);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
