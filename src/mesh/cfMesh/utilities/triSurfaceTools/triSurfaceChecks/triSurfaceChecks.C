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

#include "triSurfaceChecks.H"
#include "triSurf.H"
#include "boundBox.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "helperFunctions.H"

#ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGSurfaceChecks

# ifdef DEBUGSurfaceChecks
#include "triSurfModifier.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace triSurfaceChecks
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label checkAngles
(
    const triSurf& surf,
    labelLongList& badTriangles,
    const scalar angleTol
)
{
    badTriangles.clear();

    const scalar tol = Foam::cos(angleTol * M_PI / 180.0);

    const pointField& pts = surf.points();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(surf, tI)
    {
        const labelledTri& tri = surf[tI];

        forAll(tri, pI)
        {
            vector vn = pts[tri[tri.fcIndex(pI)]] - pts[tri[pI]];
            vn /= (mag(vn) + VSMALL);
            vector vp = pts[tri[pI]] - pts[tri[tri.rcIndex(pI)]];
            vp /= (mag(vp) + VSMALL);

            if( (vp & vn) >= tol )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                badTriangles.append(tI);

                break;
            }
        }
    }

    return badTriangles.size();
}

label checkAngles
(
    triSurf& surf,
    const word subsetName,
    const scalar angleTol
)
{
    labelLongList badTriangles;

    if( checkAngles(surf, badTriangles, angleTol) )
    {
        label setId = surf.facetSubsetIndex(subsetName);
        if( setId >= 0 )
            surf.removeFacetSubset(setId);
        setId = surf.addFacetSubset(subsetName);

        forAll(badTriangles, i)
            surf.addFacetToSubset(setId, badTriangles[i]);
    }

    return badTriangles.size();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace manifoldOps
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class surfaceNeiOp
{
    // Private data
        //- const reference to face-edges addressing
        const VRWGraph& faceEdges_;

        //- const reference to edge-faces addressing
        const VRWGraph& edgeFaces_;

public:

    // Constructors

        surfaceNeiOp
        (
            const VRWGraph& faceEdges,
            const VRWGraph& edgeFaces
        )
        :
            faceEdges_(faceEdges),
            edgeFaces_(edgeFaces)
        {}

    // Public member functions
        label size() const
        {
            return faceEdges_.size();
        }

        void operator()(const label tI, DynList<label>& neiTriangles) const
        {
            neiTriangles.clear();

            forAllRow(faceEdges_, tI, teI)
            {
                const label eI = faceEdges_(tI, teI);

                if( edgeFaces_.sizeOfRow(eI) != 2 )
                    continue;

                //- ind the neighbour triangle over the edge
                label nei = edgeFaces_(eI, 0);

                if( nei == tI )
                    nei = edgeFaces_(eI, 1);

                neiTriangles.append(nei);
            }
        }

        template<class labelListType>
        void collectGroups
        (
            std::map<label, DynList<label> >& neiGroups,
            const labelListType& elementInGroup,
            const DynList<label>& localGroupLabel
        ) const
        {

        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace manifoldOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace selectorOps
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class selectOp
{

public:

    selectOp()
    {}

    bool operator()(const label tI) const
    {
        return true;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace selectorOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


label checkSurfaceManifolds
(
    const triSurf& surf,
    labelLongList& triangleInManifold
)
{
    const VRWGraph& faceEdges = surf.facetEdges();
    const VRWGraph& edgeFaces = surf.edgeFacets();

    const label nManifolds =
        help::groupMarking
        (
            triangleInManifold,
            manifoldOps::surfaceNeiOp(faceEdges, edgeFaces),
            selectorOps::selectOp()
        );

    return nManifolds;
}

label checkSurfaceManifolds
(
    triSurf& surf,
    const word subsetPrefix
)
{
    labelLongList facetInManifold;

    const label nManifolds = checkSurfaceManifolds(surf, facetInManifold);

    if( nManifolds > 1 )
    {
        labelList groupIds(nManifolds);
        forAll(groupIds, i)
        {
            const word sName = subsetPrefix+help::scalarToText(i);
            label setId = surf.facetSubsetIndex(sName);
            if( setId >= 0 )
                surf.removeFacetSubset(setId);
            groupIds[i] = surf.addFacetSubset(sName);
        }

        forAll(facetInManifold, tI)
            surf.addFacetToSubset(groupIds[facetInManifold[tI]], tI);
    }

    return nManifolds;
}


label checkForHoles(const triSurf& surf, labelLongList& badTriangles)
{
    badTriangles.clear();

    const VRWGraph& edgeFacets = surf.edgeFacets();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(edgeFacets, eI)
    {
        if( edgeFacets.sizeOfRow(eI) == 1 )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            badTriangles.append(edgeFacets(eI, 0));
        }
    }

    return badTriangles.size();
}

label checkForHoles(triSurf& surf, const word subsetName)
{
    labelLongList trianglesNearHoles;

    if( checkForHoles(surf, trianglesNearHoles) )
    {
        label setId = surf.facetSubsetIndex(subsetName);
        if( setId >= 0 )
            surf.removeFacetSubset(setId);
        setId = surf.addFacetSubset(subsetName);

        forAll(trianglesNearHoles, i)
            surf.addFacetToSubset(setId, trianglesNearHoles[i]);
    }

    return trianglesNearHoles.size();
}

label checkForNonManifoldEdges(const triSurf& surf, labelLongList& badTriangles)
{
    badTriangles.clear();

    const VRWGraph& edgeFacets = surf.edgeFacets();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(edgeFacets, eI)
    {
        if( edgeFacets.sizeOfRow(eI) > 2 )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            {
                forAllRow(edgeFacets, eI, efI)
                    badTriangles.append(edgeFacets(eI, efI));
            }
        }
    }

    return badTriangles.size();
}

label checkForNonManifoldEdges
(
    triSurf& surf,
    const word subsetPrefix
)
{
    labelLongList trianglesNearHoles;

    if( checkForHoles(surf, trianglesNearHoles) )
    {
        label setId = surf.facetSubsetIndex(subsetPrefix);
        if( setId >= 0 )
            surf.removeFacetSubset(setId);
        setId = surf.addFacetSubset(subsetPrefix);

        forAll(trianglesNearHoles, i)
            surf.addFacetToSubset(setId, trianglesNearHoles[i]);
    }

    return trianglesNearHoles.size();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace orientationOps
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class surfaceNeiOp
{
    // Private data
        //- const reference to face-edges addressing
        const VRWGraph& faceEdges_;

        //- const reference to edge-faces addressing
        const VRWGraph& edgeFaces_;

        //- const reference to triangles
        const LongList<labelledTri>& triangles_;

public:

    // Constructors

        surfaceNeiOp
        (
            const VRWGraph& faceEdges,
            const VRWGraph& edgeFaces,
            const LongList<labelledTri>& triangles
        )
        :
            faceEdges_(faceEdges),
            edgeFaces_(edgeFaces),
            triangles_(triangles)
        {}

    // Public member functions
        label size() const
        {
            return triangles_.size();
        }

        void operator()(const label tI, DynList<label>& neiTriangles) const
        {
            neiTriangles.clear();

            const labelledTri& tri = triangles_[tI];

            forAllRow(faceEdges_, tI, teI)
            {
                const label eI = faceEdges_(tI, teI);

                if( edgeFaces_.sizeOfRow(eI) != 2 )
                    continue;

                //- ind the neighbour triangle over the edge
                label nei = edgeFaces_(eI, 0);

                if( nei == tI )
                    nei = edgeFaces_(eI, 1);

                const labelledTri& neiTri = triangles_[nei];

                //- check the orientation
                label pos(-1);
                forAll(neiTri, i)
                    if( neiTri[i] == tri[teI] )
                    {
                        pos = i;
                        break;
                    }

                if( tri[tri.fcIndex(teI)] == neiTri[neiTri.rcIndex(pos)] )
                {
                    //- triangles are of the same orientation
                    neiTriangles.append(nei);
                }
            }
        }

        template<class labelListType>
        void collectGroups
        (
            std::map<label, DynList<label> >& neiGroups,
            const labelListType& elementInGroup,
            const DynList<label>& localGroupLabel
        ) const
        {

        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace orientationOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label checkOrientation(const triSurf& surf, labelLongList& triangleInGroup)
{
    triangleInGroup.setSize(surf.size());
    triangleInGroup = -1;
    label nGroups(0);

    const VRWGraph& edgeFacets = surf.edgeFacets();
    const VRWGraph& faceEdges = surf.facetEdges();
    const LongList<labelledTri>& triangles = surf.facets();

    nGroups =
        help::groupMarking
        (
            triangleInGroup,
            orientationOps::surfaceNeiOp(faceEdges, edgeFacets, triangles),
            selectorOps::selectOp()
        );

    return nGroups;
}

label checkOrientation(triSurf& surf, const word subsetPrefix)
{
    labelLongList triangleInGroup;

    const label nGroups = checkOrientation(surf, triangleInGroup);

    if( nGroups > 1 )
    {
        labelList groupIds(nGroups);
        forAll(groupIds, i)
        {
            const word sName = subsetPrefix+help::scalarToText(i);
            label setId = surf.facetSubsetIndex(sName);
            if( setId >= 0 )
                surf.removeFacetSubset(setId);

            groupIds[i] = surf.addFacetSubset(sName);
        }

        forAll(triangleInGroup, tI)
            surf.addFacetToSubset(groupIds[triangleInGroup[tI]], tI);
    }

    return nGroups;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace connectionOps
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class surfaceNeiOp
{
    // Private data
        //- const reference to face-edges addressing
        const VRWGraph& faceEdges_;

        //- const reference to edge-faces addressing
        const VRWGraph& edgeFaces_;

public:

    // Constructors

        surfaceNeiOp
        (
            const VRWGraph& faceEdges,
            const VRWGraph& edgeFaces
        )
        :
            faceEdges_(faceEdges),
            edgeFaces_(edgeFaces)
        {}

    // Public member functions
        label size() const
        {
            return faceEdges_.size();
        }

        void operator()(const label tI, DynList<label>& neiTriangles) const
        {
            neiTriangles.clear();

            forAllRow(faceEdges_, tI, teI)
            {
                const label eI = faceEdges_(tI, teI);

                forAllRow(edgeFaces_, eI, efI)
                {
                    const label tJ = edgeFaces_(eI, efI);

                    if( tJ == tI )
                        continue;

                    neiTriangles.append(tJ);
                }
            }
        }

        template<class labelListType>
        void collectGroups
        (
            std::map<label, DynList<label> >& neiGroups,
            const labelListType& elementInGroup,
            const DynList<label>& localGroupLabel
        ) const
        {

        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace orientationOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label checkDisconnectedParts
(
    const triSurf& surf,
    labelLongList& triangleInRegion
)
{
    triangleInRegion.setSize(surf.size());
    triangleInRegion = -1;
    label nGroups(0);

    const VRWGraph& edgeFacets = surf.edgeFacets();
    const VRWGraph& faceEdges = surf.facetEdges();

    nGroups =
        help::groupMarking
        (
            triangleInRegion,
            connectionOps::surfaceNeiOp(faceEdges, edgeFacets),
            selectorOps::selectOp()
        );

    return nGroups;
}

label checkDisconnectedParts(triSurf& surf, const word subsetPrefix)
{
    labelLongList triangleInRegion;

    const label nGroups = checkDisconnectedParts(surf, triangleInRegion);

    if( nGroups > 1 )
    {
        labelList groupIds(nGroups);
        forAll(groupIds, i)
        {
            const word sName = subsetPrefix+help::scalarToText(i);
            label setId = surf.facetSubsetIndex(sName);
            if( setId >= 0 )
                surf.removeFacetSubset(setId);

            groupIds[i] = surf.addFacetSubset(sName);
        }

        forAll(triangleInRegion, tI)
            surf.addFacetToSubset(groupIds[triangleInRegion[tI]], tI);
    }

    return nGroups;
}

void calculateBoundingBox(const triSurf& surf, boundBox& bb)
{
    bb.min() = Foam::min(surf.points());
    bb.max() = Foam::max(surf.points());
}

label checkCollocatedPoints
(
    const triSurf& surf,
    labelLongList& collocatedPoints,
    const scalar distTol
)
{
    collocatedPoints.clear();

    meshOctree octree(surf);
    meshOctreeCreator(octree).createOctreeWithRefinedBoundary(20, 30);

    const pointField& pts = surf.points();

    boolList collocated(pts.size(), false);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(collocated, pI)
    {
        const point& p = pts[pI];

        boundBox bb(p, p);

        bb.min() -= point(distTol, distTol, distTol);
        bb.max() += point(distTol, distTol, distTol);

        DynList<label> leavesInBox;
        octree.findLeavesContainedInBox(bb, leavesInBox);

        forAll(leavesInBox, i)
        {
            const label leafI = leavesInBox[i];

            DynList<label> trianglesInBox;

            octree.containedTriangles(leafI, trianglesInBox);

            forAll(trianglesInBox, j)
            {
                const label triJ = trianglesInBox[j];
                const labelledTri& nt = surf[triJ];

                forAll(nt, tpI)
                {
                    if( nt[tpI] == pI )
                        continue;

                    if( magSqr(pts[nt[tpI]] - p) < sqr(distTol) )
                    {
                        collocated[pI] = true;
                        collocated[nt[tpI]] = true;
                    }
                }
            }
        }
    }

    forAll(collocated, pI)
        if( collocated[pI] )
            collocatedPoints.append(pI);

    return collocatedPoints.size();
}

label checkCollocatedPoints
(
    triSurf& surf,
    const word subsetName,
    const scalar distTol
)
{
    labelLongList collocatedPoints;

    if( checkCollocatedPoints(surf, collocatedPoints, distTol) )
    {
        label setId = surf.pointSubsetIndex(subsetName);
        if( setId >= 0 )
            surf.removePointSubset(setId);
        setId = surf.addPointSubset(subsetName);

        forAll(collocatedPoints, i)
            surf.addPointToSubset(setId, collocatedPoints[i]);
    }

    return collocatedPoints.size();
}

label checkSelfIntersections
(
    const triSurf& surf,
    labelLongList& badFaces,
    const scalar tol
)
{
    badFaces.clear();

    meshOctree octree(surf);
    meshOctreeCreator(octree).createOctreeWithRefinedBoundary(20, 30);

    const pointField& pts = surf.points();

    boolList intersected(surf.size(), false);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(intersected, tI)
    {
        const labelledTri& tri = surf[tI];

        const triangle<point, point> currTri
        (
            pts[tri[0]],
            pts[tri[1]],
            pts[tri[2]]
        );

        boundBox bb(pts[tri[0]], pts[tri[0]]);
        for(label i=1;i<3;++i)
        {
            bb.min() = Foam::min(bb.min(), pts[tri[i]]);
            bb.max() = Foam::max(bb.max(), pts[tri[i]]);
        }

        bb.min() -= point(tol, tol, tol);
        bb.max() += point(tol, tol, tol);

        DynList<label> leavesInBox;
        octree.findLeavesContainedInBox(bb, leavesInBox);

        forAll(leavesInBox, i)
        {
            const label leafI = leavesInBox[i];

            DynList<label> trianglesInBox;

            octree.containedTriangles(leafI, trianglesInBox);

            forAll(trianglesInBox, j)
            {
                const label triJ = trianglesInBox[j];
                const labelledTri& nt = surf[triJ];

                if( tI >= triJ )
                    continue;

                DynList<label, 3> sharedPoints;
                forAll(nt, pJ)
                {
                    forAll(tri, pI)
                    if( tri[pI] == nt[pJ] )
                        sharedPoints.append(pJ);
                }

                if( sharedPoints.size() >= 2 )
                {
                    //- triangles share an edge and cannot self-intersect
                    continue;
                }
                else if( sharedPoints.size() == 1 )
                {
                    //- check if the opposite edge intersects the triangle
                    const point& s = pts[nt[(sharedPoints[0]+1)%3]];
                    const point& e = pts[nt[(sharedPoints[0]+2)%3]];

                    point intersection;
                    const bool lineIntersection =
                        help::triLineIntersection(currTri, s, e, intersection);

                    if( lineIntersection )
                    {
                        intersected[tI] = true;
                        intersected[triJ] = true;
                    }
                }
                else
                {
                    //- triangles do not share any vertices
                    const triangle<point, point> neiTri
                    (
                        pts[nt[0]],
                        pts[nt[1]],
                        pts[nt[2]]
                    );

                    const bool intersect =
                        help::doTrianglesIntersect
                        (
                            currTri,
                            neiTri,
                            tol
                        );

                    if( intersect )
                    {
                        intersected[tI] = true;
                        intersected[triJ] = true;

                        # ifdef DEBUGSurfaceChecks
                        const label sId =
                            const_cast<triSurf&>(surf).addFacetSubset
                            (
                                "noCommonIntersect_"+help::scalarToText(tI)+
                                "_"+help::scalarToText(triJ)
                            );
                        const_cast<triSurf&>(surf).addFacetToSubset(sId, tI);
                        const_cast<triSurf&>(surf).addFacetToSubset(sId, triJ);
                        # endif
                    }
                }
            }
        }
    }

    forAll(intersected, tI)
    {
        if( intersected[tI] )
            badFaces.append(tI);
    }

    return badFaces.size();
}

label checkSelfIntersections
(
    triSurf& surf,
    const word subsetName,
    const scalar tol
)
{
    labelLongList badFacets;

    if( checkSelfIntersections(surf, badFacets, tol) )
    {
        label setId = surf.facetSubsetIndex(subsetName);
        if( setId >= 0 )
            surf.removeFacetSubset(setId);
        setId = surf.addFacetSubset(subsetName);

        forAll(badFacets, i)
            surf.addFacetToSubset(setId, badFacets[i]);
    }

    return badFacets.size();
}

label checkOverlaps
(
    const triSurf& surf,
    labelLongList& badFaces,
    const scalar tol,
    const scalar angleTol
)
{
    badFaces.clear();

    meshOctree octree(surf);
    meshOctreeCreator(octree).createOctreeWithRefinedBoundary(20, 30);

    const scalar cosVal = Foam::cos(angleTol * M_PI / 180.0);

    const pointField& pts = surf.points();

    boolList intersected(surf.size(), false);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(intersected, tI)
    {
        const labelledTri& tri = surf[tI];

        const triangle<point, point> currTri
        (
            pts[tri[0]],
            pts[tri[1]],
            pts[tri[2]]
        );

        boundBox bb(pts[tri[0]], pts[tri[0]]);
        for(label i=1;i<3;++i)
        {
            bb.min() = Foam::min(bb.min(), pts[tri[i]]);
            bb.max() = Foam::max(bb.max(), pts[tri[i]]);
        }

        bb.min() -= point(tol, tol, tol);
        bb.max() += point(tol, tol, tol);

        DynList<label> leavesInBox;
        octree.findLeavesContainedInBox(bb, leavesInBox);

        forAll(leavesInBox, i)
        {
            const label leafI = leavesInBox[i];

            DynList<label> trianglesInBox;

            octree.containedTriangles(leafI, trianglesInBox);

            forAll(trianglesInBox, j)
            {
                const label triJ = trianglesInBox[j];
                const labelledTri& nt = surf[triJ];

                if( tI >= triJ )
                    continue;

                const triangle<point, point> neiTri
                (
                    pts[nt[0]],
                    pts[nt[1]],
                    pts[nt[2]]
                );

                DynList<point> commonPolygon;

                const bool intersect =
                    help::doTrianglesOverlap
                    (
                        currTri,
                        neiTri,
                        commonPolygon,
                        tol,
                        cosVal
                    );

                if( intersect )
                {
                    intersected[tI] = true;
                    intersected[triJ] = true;

                    # ifdef DEBUGSurfaceChecks
                    const label sId =
                        const_cast<triSurf&>(surf).addFacetSubset
                        (
                            "noCommonOverlap_"+help::scalarToText(tI)+
                            "_"+help::scalarToText(triJ)
                        );
                    const_cast<triSurf&>(surf).addFacetToSubset(sId, tI);
                    const_cast<triSurf&>(surf).addFacetToSubset(sId, triJ);

                    triSurf newSurf;
                    geometricSurfacePatchList& patches =
                        triSurfModifier(newSurf).patchesAccess();
                    patches.setSize(1);
                    patches[0].name() = "patch0";
                    forAll(commonPolygon, i)
                        newSurf.appendVertex(commonPolygon[i]);
                    for(label i=0;i<commonPolygon.size()-2;++i)
                        newSurf.appendTriangle(labelledTri(0, i+1, i+2, 0));
                    newSurf.writeSurface
                    (
                        "overlap_"+help::scalarToText(tI)+
                        "_"+help::scalarToText(triJ)+".fms"
                    );
                    # endif
                }
            }
        }
    }

    forAll(intersected, tI)
    {
        if( intersected[tI] )
            badFaces.append(tI);
    }

    return badFaces.size();
}

label checkOverlaps
(
    triSurf& surf,
    const word subsetName,
    const scalar tol,
    const scalar angleTol
)
{
    labelLongList badFaces;

    if( checkOverlaps(surf, badFaces, tol, angleTol) )
    {
        label setId = surf.facetSubsetIndex(subsetName);
        if( setId >= 0 )
            surf.removeFacetSubset(setId);
        setId = surf.addFacetSubset(subsetName);

        forAll(badFaces, i)
            surf.addFacetToSubset(setId, badFaces[i]);
    }

    return badFaces.size();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace triSurfaceChecks

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
