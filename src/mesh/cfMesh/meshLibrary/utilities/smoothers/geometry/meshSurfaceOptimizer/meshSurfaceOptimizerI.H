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

#include "demandDrivenData.H"
#include "meshSurfaceOptimizer.H"
#include "plane.H"
#include "Map.H"
#include "surfaceOptimizer.H"
#include "helperFunctions.H"
#include "partTriMeshSimplex.H"

//#define DEBUGSmooth

# ifdef DEBUGSmooth
#include "triSurf.H"
#include "triSurfModifier.H"
#include "helperFunctions.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const partTriMesh& meshSurfaceOptimizer::triMesh() const
{
    if( !triMeshPtr_ )
        calculateTrianglesAndAddressing();

    return *triMeshPtr_;
}

inline void meshSurfaceOptimizer::updateTriMesh(const labelLongList& selPoints)
{
    if( !triMeshPtr_ )
        FatalErrorIn
        (
            "inline void meshSurfaceOptimizer::updateTriMesh"
            "(const labelLongList&)"
        ) << "triMeshPtr_ is not allocated " << abort(FatalError);

    triMeshPtr_->updateVertices(selPoints);
}

inline void meshSurfaceOptimizer::updateTriMesh()
{
    if( !triMeshPtr_ )
        FatalErrorIn
        (
            "inline void meshSurfaceOptimizer::updateTriMesh()"
        ) << "triMeshPtr_ is not allocated " << abort(FatalError);

    triMeshPtr_->updateVertices();
}

inline bool meshSurfaceOptimizer::transformIntoPlane
(
    const label bpI,
    const plane& pl,
    vector& vecX,
    vector& vecY,
    DynList<point>& pts,
    DynList<triFace>& trias
) const
{
    # ifdef DEBUGSmooth
    Info << "Transforming boundary node " << bpI << endl;
    # endif

    const partTriMesh& triMesh = this->triMesh();
    const label triPointI = triMesh.meshSurfacePointLabelInTriMesh()[bpI];
    partTriMeshSimplex simplex(triMesh, triPointI);

    const DynList<point, 32>& sPts = simplex.pts();

    //- create vecX and vecY
    const point& p = simplex.centrePoint();
    pts.setSize(0);
    bool found(false);
    forAll(sPts, pI)
    {
        const point sp = pl.nearestPoint(sPts[pI]);
        const scalar d = mag(sp - p);
        if( d > VSMALL )
        {
            vecX = sp - p;
            vecX /= d;
            vecY = pl.normal() ^ vecX;
            vecY /= mag(vecY);
            found = true;
            break;
        }
    }

    if( !found )
        return false;

    trias = simplex.triangles();

    # ifdef DEBUGSmooth
    Info << "VecX " << vecX << endl;
    Info << "vecY " << vecY << endl;
    # endif

    //- transform the vertices
    pts.setSize(sPts.size());

    forAll(sPts, pI)
    {
        const point pt
        (
            ((sPts[pI] - p) & vecX),
            ((sPts[pI] - p) & vecY),
            0.0
        );

        pts[pI] = pt;
    }

    # ifdef DEBUGSmooth
    Info << "Original triangles " << endl;
    forAll(simplex.triangles(), triI)
        Info << "Tri " << triI << " is " << simplex.triangles()[triI] << endl;
    Info << "Transformed triangles are " << trias << endl;
    Info << "Transformed vertices " << pts << endl;

    triSurf surf;
    triSurfModifier sMod(surf);
    pointField& sPoints = sMod.pointsAccess();
    sPoints.setSize(pts.size());
    forAll(sPoints, i)
        sPoints[i] = pts[i];
    LongList<labelledTri>& sTrias = sMod.facetsAccess();
    sTrias.setSize(trias.size());
    forAll(sTrias, i)
    {
        labelledTri& tf = sTrias[i];
        tf[0] = trias[i][0];
        tf[1] = trias[i][1];
        tf[2] = trias[i][2];

        tf.region() = 0;
    }
    sMod.patchesAccess().setSize(1);
    sMod.patchesAccess()[0].name() = "bnd";
    sMod.patchesAccess()[0].geometricType() = "patch";

    fileName sName("bndSimplex_");
    sName += help::scalarToText(bpI);
    sName += ".stl";
    surf.writeSurface(sName);
    # endif

    return found;
}

inline bool meshSurfaceOptimizer::transformIntoPlaneParallel
(
    const label bpI,
    const plane& pl,
    const std::map<label, DynList<parTriFace> >& m,
    vector& vecX,
    vector& vecY,
    DynList<point>& pts,
    DynList<triFace>& trias
) const
{
    /*
    # ifdef DEBUGSmooth
    Info << "Transforming boundary node at parallel interface " << bpI << endl;
    # endif

    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    const labelList& globalPointLabel =
        surfaceEngine_.globalBoundaryPointLabel();

    std::map<label, DynList<parTriFace> >::const_iterator
        pIter = m.find(globalPointLabel[bpI]);

    if( pIter == m.end() )
    {
        FatalErrorIn
        (
            "(const label, const plane&,"
            " const std::map<label, DynList<parTriFace> >&,"
            " vector&, vector&, DynList<point>&, DynList<triFace>&) const"
        ) << " Cannot find node " << globalPointLabel[bpI]
          << " in the map" << abort(FatalError);
    }

    const DynList<parTriFace>& pFcs = pIter->second;

    # ifdef DEBUGSmooth
    if( globalPointLabel[bpI] != pFcs[0].globalLabelOfPoint(0) )
        FatalError << "Wrong points supplied!" << abort(FatalError);
    Info << "Triangles containing this node are " << pFcs << endl;
    # endif

    //- create vecX and vecY
    const point& p = points[bPoints[bpI]];
    pts.setSize(0);
    bool found(false);

    forAll(pFcs, tI)
    {
        const point sp = pl.nearestPoint(pFcs[tI].trianglePoints().centre());
        const scalar d = mag(sp - p);
        if( d > VSMALL )
        {
            vecX = sp - p;
            vecX /= d;
            vecY = pl.normal() ^ vecX;
            vecY /= mag(vecY);
            found = true;
            break;
        }
    }

    if( !found )
        return false;

    # ifdef DEBUGSmooth
    Info << "VecX " << vecX << endl;
    Info << "vecY " << vecY << endl;
    # endif

    //- transform the vertices
    //label counter(0);
    Map<label> pointMapping;
    trias.setSize(pFcs.size());

    //- tranfer triangle from other proc first
    forAll(pFcs, triI)
    {
        const parTriFace& tri = pFcs[triI];

        triFace& tf = trias[triI];

        for(label pI=0;pI<3;++pI)
        {
            const Map<label>::iterator fnd =
                pointMapping.find(tri.globalLabelOfPoint(pI));

            if( fnd != pointMapping.end() )
            {
                tf[pI] = fnd();
            }
            else
            {
                point tp;
                if( pI == 0 )
                {
                    tp = tri.trianglePoints().a();
                }
                else if( pI == 1 )
                {
                    tp = tri.trianglePoints().b();
                }
                else
                {
                    tp = tri.trianglePoints().c();
                }

                const point pt
                (
                    ((tp - p) & vecX),
                    ((tp - p) & vecY),
                    0.0
                );

                tf[pI] = pts.size();
                pointMapping.insert(tri.globalLabelOfPoint(pI), pts.size());
                pts.append(pt);
            }
        }
    }

    return found;
    */

    return false;
}

inline point meshSurfaceOptimizer::newPositionLaplacian
(
    const label bpI,
    const bool transformIntoPlane
) const
{
    const VRWGraph& pPoints = surfaceEngine_.pointPoints();
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    vector newP(vector::zero);
    if( transformIntoPlane )
    {
        const vector& pNormal = surfaceEngine_.pointNormals()[bpI];

        if( magSqr(pNormal) < VSMALL )
            return points[bPoints[bpI]];

        plane pl(points[bPoints[bpI]], pNormal);

        DynList<point> projectedPoints;
        projectedPoints.setSize(pPoints.sizeOfRow(bpI));
        forAllRow(pPoints, bpI, pI)
            projectedPoints[pI] =
                pl.nearestPoint(points[bPoints[pPoints(bpI, pI)]]);

        forAll(projectedPoints, pI)
            newP += projectedPoints[pI];

        newP /= projectedPoints.size();
    }
    else
    {
        forAllRow(pPoints, bpI, pI)
            newP += points[bPoints[pPoints(bpI, pI)]];

        newP /= pPoints.sizeOfRow(bpI);
    }

    return newP;
}

inline point meshSurfaceOptimizer::newPositionLaplacianFC
(
    const label bpI,
    const bool transformIntoPlane
) const
{
    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();
    const pointFieldPMG& points = surfaceEngine_.points();
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    vector newP(vector::zero);

    if( transformIntoPlane )
    {
        const vector& pNormal = surfaceEngine_.pointNormals()[bpI];

        if( magSqr(pNormal) < VSMALL )
            return points[bPoints[bpI]];

        plane pl(points[bPoints[bpI]], pNormal);

        DynList<point> projectedPoints;
        projectedPoints.setSize(pointFaces.sizeOfRow(bpI));
        forAllRow(pointFaces, bpI, pfI)
            projectedPoints[pfI] =
                pl.nearestPoint(faceCentres[pointFaces(bpI, pfI)]);

        forAll(projectedPoints, pI)
            newP += projectedPoints[pI];

        newP /= projectedPoints.size();
    }
    else
    {
        forAllRow(pointFaces, bpI, pfI)
            newP += faceCentres[pointFaces(bpI, pfI)];

        newP /= pointFaces.sizeOfRow(bpI);
    }

    return newP;
}

inline point meshSurfaceOptimizer::newPositionLaplacianWFC
(
    const label bpI,
    const bool transformIntoPlane
) const
{
    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();
    const pointFieldPMG& points = surfaceEngine_.points();
    const vectorField& faceAreas = surfaceEngine_.faceNormals();
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    vector newP(vector::zero);

    if( transformIntoPlane )
    {
        const vector& pNormal = surfaceEngine_.pointNormals()[bpI];

        if( magSqr(pNormal) < VSMALL )
            return points[bPoints[bpI]];

        plane pl(points[bPoints[bpI]], pNormal);

        DynList<point> projectedPoints;
        projectedPoints.setSize(pointFaces.sizeOfRow(bpI));
        forAllRow(pointFaces, bpI, pfI)
            projectedPoints[pfI] =
                pl.nearestPoint(faceCentres[pointFaces(bpI, pfI)]);

        scalar sumW(0.0);
        forAll(projectedPoints, pI)
        {
            const label bfI = pointFaces(bpI, pI);
            const scalar w = (faceAreas[bfI] & faceAreas[bfI]);
            newP += w * projectedPoints[pI];
            sumW += w;
        }

        newP /= Foam::max(sumW, VSMALL);
    }
    else
    {
        scalar sumW(0.0);
        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            const scalar w = (faceAreas[bfI] & faceAreas[bfI]);
            newP += w * faceCentres[pointFaces(bpI, pfI)];
            sumW += w;
        }

        newP /= Foam::max(sumW, VSMALL);
    }

    return newP;
}

inline point meshSurfaceOptimizer::newPositionSurfaceOptimizer
(
    const label bpI,
    const scalar tol
) const
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    # ifdef DEBUGSmooth
    Info << "Smoothing boundary node " << bpI << endl;
    Info << "Node label in the mesh is " << bPoints[bpI] << endl;
    Info << "Point coordinates " << points[bPoints[bpI]] << endl;
    # endif

    //- project vertices onto the plane
    const vector& pNormal = surfaceEngine_.pointNormals()[bpI];
    if( magSqr(pNormal) < VSMALL )
        return points[bPoints[bpI]];

    const plane pl(points[bPoints[bpI]], pNormal);

    DynList<point> pts;
    DynList<triFace> trias;
    vector vecX, vecY;
    bool success = this->transformIntoPlane(bpI, pl, vecX, vecY, pts, trias);
    if( !success )
    {
        Warning << "Cannot transform to plane" << endl;
        return points[bPoints[bpI]];
    }

    surfaceOptimizer so(pts, trias);
    const point newPoint = so.optimizePoint(tol);

    const point newP
    (
        points[bPoints[bpI]] +
        vecX * newPoint.x() +
        vecY * newPoint.y()
    );

    if( help::isnan(newP) || help::isinf(newP) )
    {
        WarningIn
        (
            "inline point meshSurfaceOptimizer::newPositionSurfaceOptimizer"
            "(const label, const scalar) const"
        ) << "Cannot move point " << bpI << endl;

        return points[bPoints[bpI]];
    }

    return newP;
}

inline point meshSurfaceOptimizer::newEdgePositionLaplacian
(
    const label bpI
) const
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const edgeList& edges = surfaceEngine_.edges();
    const VRWGraph& bpEdges = surfaceEngine_.boundaryPointEdges();
    const pointFieldPMG& points = surfaceEngine_.points();

    const labelHashSet& featureEdges = partitionerPtr_->featureEdges();

    DynList<label> edgePoints;

    forAllRow(bpEdges, bpI, i)
    {
        const label beI = bpEdges(bpI, i);

        if( featureEdges.found(beI) )
        {
            edgePoints.append(edges[beI].otherVertex(bPoints[bpI]));
        }
    }

    if( edgePoints.size() != 2 )
        return points[bPoints[bpI]];

    # ifdef DEBUGSearch
    Info << "Edge points " << edgePoints << endl;
    # endif

    vector pos(vector::zero);
    forAll(edgePoints, epI)
        pos += points[edgePoints[epI]];

    pos /= edgePoints.size();

    return pos;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
