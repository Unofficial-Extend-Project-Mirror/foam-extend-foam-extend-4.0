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

#include "error.H"
#include "helperFunctionsGeometryQueries.H"
#include "helperFunctionsTopologyManipulation.H"
#include "edgeList.H"
#include "pointField.H"
#include "boolList.H"
#include "triSurf.H"
#include "matrix3D.H"
#include "HashSet.H"
#include "tetrahedron.H"
#include "boundBox.H"
#include "Map.H"

//#define DEBUG_pMesh

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

namespace help
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

template<class ListType>
inline bool isnan(const ListType& l)
{
    forAll(l, i)
        if( l[i] != l[i] )
            return true;

    return false;
}

template<class ListType>
bool isinf(const ListType& l)
{
    forAll(l, i)
        if( (l[i] < -VGREAT) || (l[i] > VGREAT) )
            return true;

    return false;
}

template<class Face1, class Face2>
inline bool isSharedEdgeConvex
(
    const pointField& points,
    const Face1& f1,
    const Face2& f2
)
{
    DynList<label, 3> triOwn(3);
    DynList<label, 3> triNei(3);

    forAll(f1, pI)
    {
        label pos(-1);
        forAll(f2, pJ)
            if( f2[pJ] == f1[pI] )
            {
                pos = pJ;
                break;
            }

        if( pos < 0 )
            continue;

        triNei[0] = f2[pos];
        triNei[1] = f2[f2.fcIndex(pos)];
        triNei[2] = f2[f2.rcIndex(pos)];

        triOwn[0] = f1[pI];
        triOwn[1] = f1[f1.fcIndex(pI)];
        triOwn[2] = f1[f1.rcIndex(pI)];

        scalar vol(0.0);

        forAll(triOwn, pJ)
        {
            if( !triNei.contains(triOwn[pJ]) )
            {
                tetrahedron<point, point> tet
                (
                    points[triNei[0]],
                    points[triNei[1]],
                    points[triNei[2]],
                    points[triOwn[pJ]]
                );

                vol = tet.mag();
                break;
            }
        }

        if( vol > -VSMALL )
            return false;
    }

    return true;
}

template<class Face1, class Face2>
inline scalar angleBetweenFaces
(
    const pointField& points,
    const Face1& f1,
    const Face2& f2
)
{
    DynList<label, 3> triOwn(3);
    DynList<label, 3> triNei(3);

    scalar angle(0.0);
    label counter(0);

    forAll(f1, pI)
    {
        label pos(-1);
        forAll(f2, pJ)
            if( f2[pJ] == f1[pI] )
            {
                pos = pJ;
                break;
            }

        if( pos < 0 )
            continue;

        triNei[0] = f2[pos];
        triNei[1] = f2[f2.fcIndex(pos)];
        triNei[2] = f2[f2.rcIndex(pos)];

        triOwn[0] = f1[pI];
        triOwn[1] = f1[f1.fcIndex(pI)];
        triOwn[2] = f1[f1.rcIndex(pI)];

        scalar vol(0.0);

        forAll(triOwn, pJ)
        {
            if( !triNei.contains(triOwn[pJ]) )
            {
                tetrahedron<point, point> tet
                (
                    points[triNei[0]],
                    points[triNei[1]],
                    points[triNei[2]],
                    points[triOwn[pJ]]
                );

                vol = tet.mag();
                break;
            }
        }

        vector nOwn
        (
            (points[triOwn[1]] - points[triOwn[0]]) ^
            (points[triOwn[2]] - points[triOwn[0]])
        );
        nOwn /= (mag(nOwn) + VSMALL);

        vector nNei
        (
            (points[triNei[1]] - points[triNei[0]]) ^
            (points[triNei[2]] - points[triNei[0]])
        );
        nNei /= (mag(nNei) + VSMALL);

        const scalar dot = Foam::max(-1.0, Foam::min(nOwn & nNei, 1.0));

        if( vol > -VSMALL )
        {
            //- the angle is in the interval [Pi, 2Pi>
            const scalar ang = Foam::acos(dot);

            angle += ang + M_PI;
            ++counter;
        }
        else
        {
            //- the angle is in the interval [0, Pi>
            const scalar ang = Foam::acos(-dot);

            angle += ang;
            ++counter;
        }
    }

    if( counter == 0 )
    {
        FatalErrorIn
        (
            "scalar angleBetweenFaces"
            "(const pointField&, const face&, const face&)"
        ) << "Faces " << f1 << " and " << f2
          << " do no share an edge" << abort(FatalError);
    }

    return angle / counter;
}

inline faceList mergePatchFaces
(
    const List< DynList<label> >& pfcs,
    const pointField& polyPoints
)
{
    //- merge faces which share a common edge
    faceList patchFaces(pfcs.size());
    label counter(0);
    forAll(pfcs, faceI)
        if( pfcs[faceI].size() > 2 )
        {
            const DynList<label>& f = pfcs[faceI];
            face f_(f.size());
            forAll(f_, fJ)
                f_[fJ] = f[fJ];

            patchFaces[counter++] = f_;
        }

    patchFaces.setSize(counter);

    bool merged;
    do
    {
        faceList mergedFaces(patchFaces.size());
        boolList currentlyMerged(patchFaces.size(), false);

        counter = 0;
        merged = false;

        for(label nI=0;nI<(patchFaces.size()-1);nI++)
        {
            vector n0 = patchFaces[nI].normal(polyPoints);
            n0 /= mag(n0);

            for(label nJ=nI+1;nJ<patchFaces.size();nJ++)
            {
                vector n1 = patchFaces[nI].normal(polyPoints);
                n1 /= mag(n1);
                if(
                    help::shareAnEdge(patchFaces[nI], patchFaces[nJ]) &&
                    mag(n0 & n1) > 0.95
                )
                {
                    merged = true;
                    currentlyMerged[nI] = currentlyMerged[nJ] = true;
                    mergedFaces[counter++] =
                        help::mergeTwoFaces
                        (
                            patchFaces[nI],
                            patchFaces[nJ]
                        );

                    break;
                }
            }

            if( merged ) break;
        }

        forAll(patchFaces, pfI)
            if( !currentlyMerged[pfI] )
                mergedFaces.newElmt(counter++) = patchFaces[pfI];

        if( merged )
        {
            patchFaces.setSize(counter);
            for(label k=0;k<counter;k++)
                patchFaces[k] = mergedFaces[k];
        }

    } while( merged );

    return patchFaces;
}

inline bool vertexOnLine(const point& p, const edge& e, const pointField& ep)
{
    vector v = e.vec(ep);
    v /= mag(v);

    vector pv = p - ep[e.start()];
    pv /= mag(pv);

    if( mag(pv & v) > (1.0-SMALL) )
        return true;

    return false;
}

inline bool vertexInPlane(const point& p, const plane& pl)
{
    const vector& n = pl.normal();
    const point& fp = pl.refPoint();

    vector d = p - fp;
    if( mag(d) > VSMALL )
        d /= mag(d);

    if( mag(d & n) < SMALL )
        return true;

    return false;
}

inline bool planeIntersectsEdge
(
    const point& start,
    const point& end,
    const plane& pl,
    point& intersection
)
{
    const vector v = end - start;

    const vector& n = pl.normal();
    const point& fp = pl.refPoint();

    if( (n & (v/mag(v))) < SMALL )
        return false;

    const scalar t((n & (fp - start)) / (n & v));

    if( (t > -SMALL) && (t < (1.0+SMALL)) )
    {
        intersection = start + v * t;
        return true;
    }

    return false;
}

inline bool pointInTetrahedron
(
    const point& p,
    const tetrahedron<point, point>& tet
)
{
    const vector v0 = tet.a() - tet.d();
    const vector v1 = tet.b() - tet.d();
    const vector v2 = tet.c() - tet.d();
    const vector sp = p - tet.d();

    matrix3D mat;
    FixedList<scalar, 3> source;
    for(label i=0;i<3;++i)
    {
        mat[i][0] = v0[i];
        mat[i][1] = v1[i];
        mat[i][2] = v2[i];
        source[i] = sp[i];
    }

    //- check the determinant of the transformation
    const scalar det = mat.determinant();

    if( mag(det) < VSMALL )
        return false;

    //- get the coordinates of the point in the barycentric corrdinate system
    const scalar u0 = mat.solveFirst(source);

    if( (u0 < -SMALL) || (u0 > (1.0+SMALL)) )
        return false;

    const scalar u1 = mat.solveSecond(source);

    if( (u1 < -SMALL) || ((u0+u1) > (1.0+SMALL)) )
        return false;

    const scalar u2 = mat.solveThird(source);

    if( (u2 < -SMALL) || (u2 > (1.0+SMALL)) )
        return false;


    const scalar u3 = 1.0 - u0 - u1 - u2;

    if( (u3 < -SMALL) || (u3 > (1.0+SMALL)) )
        return false;

    return true;
}

inline bool nearestEdgePointToTheLine
(
    const point& edgePoint0,
    const point& edgePoint1,
    const point& lp0,
    const point& lp1,
    point& nearestOnEdge,
    point& nearestOnLine
)
{
    const vector v = lp1 - lp0;
    const vector d = lp0 - edgePoint0;
    const vector e = edgePoint1 - edgePoint0;

    const scalar vMag = mag(v);
    if( vMag < VSMALL )
        return false;

    const scalar eMag = mag(e);
    if( eMag < VSMALL )
    {
        nearestOnEdge = edgePoint0;
        nearestOnLine = nearestPointOnTheEdge(lp0, lp1, nearestOnEdge);
        return true;
    }

    if( mag((v/vMag) & (e/eMag)) > (1.0 - SMALL) )
        return false;

    tensor mat(tensor::zero);
    mat.xx() = (v&v);
    mat.xy() = mat.yx() = -1.0 * (v&e);
    mat.yy() = (e&e);
    mat.zz() = SMALL;

    vector source(vector::zero);
    source[0] = -1.0 * (d&v);
    source[1] = (d&e);

    const vector sol = (inv(mat) & source);

    nearestOnLine = lp0 + v * sol[0];
    if( sol[1] > 1.0 )
    {
        nearestOnEdge = edgePoint1;
    }
    else if( sol[1] < 0.0 )
    {
        nearestOnEdge = edgePoint0;
    }
    else
    {
        nearestOnEdge = edgePoint0 + e * sol[1];
    }

    return true;
}

inline point nearestPointOnTheTriangle
(
    const label tI,
    const triSurf& surface,
    const point& p
)
{
    const labelledTri& ltri = surface[tI];
    const pointField& points = surface.points();

    const vector v0 = points[ltri[1]] - points[ltri[0]];
    const vector v1 = points[ltri[2]] - points[ltri[0]];
    const vector v2 = p - points[ltri[0]];

    const scalar dot00 = (v0 & v0);
    const scalar dot01 = (v0 & v1);
    const scalar dot02 = (v0 & v2);
    const scalar dot11 = (v1 & v1);
    const scalar dot12 = (v1 & v2);

    // Compute barycentric coordinates
    const scalar det = dot00 * dot11 - dot01 * dot01;

    if( mag(det) < VSMALL )
    {
        point nearest(p);
        scalar dist(VGREAT);
        const edgeList edges = ltri.edges();
        forAll(edges, eI)
        {
            const edge& e = edges[eI];
            const point np =
                nearestPointOnTheEdge
                (
                    points[e.start()],
                    points[e.end()],
                    p
                );

            if( magSqr(p - np) < dist )
            {
                nearest = np;
                dist = magSqr(p - np);
            }
        }

        return nearest;
    }

    const scalar u = (dot11 * dot02 - dot01 * dot12) / det;
    const scalar v = (dot00 * dot12 - dot01 * dot02) / det;

    const point pProj = points[ltri[0]] + u * v0 + v * v1;

    //- Check if point is in triangle
    if( (u >= -SMALL) && (v >= -SMALL) && ((u + v) <= (1.0+SMALL)) )
    {
        return pProj;
    }
    else
    {
        if( u < -SMALL )
        {
            const edge e(ltri[0], ltri[2]);
            const scalar ed = ((pProj - points[e.start()]) & v1) / magSqr(v1);

            if( ed > 1.0 )
            {
                return points[e.end()];
            }
            else if( ed < 0.0 )
            {
                return points[e.start()];
            }
            else
            {
                return (points[e.start()] + v1 * ed);
            }
        }
        else if( v < -SMALL )
        {
            const edge e(ltri[0], ltri[1]);
            const scalar ed = ((pProj - points[e.start()]) & v0) / magSqr(v0);

            if( ed > 1.0 )
            {
                return points[e.end()];
            }
            else if( ed < 0.0 )
            {
                return points[e.start()];
            }
            else
            {
                return (points[e.start()] + v0 * ed);
            }
        }
        else
        {
            const edge e(ltri[2], ltri[1]);
            const vector ev = e.vec(points);
            const scalar ed = ((pProj - points[e.start()]) & ev) / magSqr(ev);

            if( ed > 1.0 )
            {
                return points[e.end()];
            }
            else if( ed < 0.0 )
            {
                return points[e.start()];
            }
            else
            {
                return (points[e.start()] + ev * ed);
            }
        }
    }

    return pProj;
}

inline bool triLineIntersection
(
    const triangle<point, point>& tria,
    const point& lineStart,
    const point& lineEnd,
    point& intersection
)
{
    const point& p0 = tria.a();
    const vector v(lineStart - lineEnd);
    const vector v0 = tria.b() - p0;
    const vector v1 = tria.c() - p0;
    const vector sp = lineStart - p0;

    matrix3D mat;
    FixedList<scalar, 3> source;
    for(label i=0;i<3;++i)
    {
        mat[i][0] = v0[i];
        mat[i][1] = v1[i];
        mat[i][2] = v[i];
        source[i] = sp[i];
    }

    const scalar det = mat.determinant();

    if( mag(det) < SMALL )
        return false;

    const scalar t = mat.solveThird(source);

    if( (t < -SMALL) || (t > (1.0+SMALL)) )
        return false;

    const scalar u0 = mat.solveFirst(source);

    if( u0 < -SMALL )
        return false;

    const scalar u1 = mat.solveSecond(source);

    if( (u1 < -SMALL) || ((u0+u1) > (1.0+SMALL)) )
        return false;

    intersection = lineStart - t * v;
    return true;
}

inline bool triLineIntersection
(
    const triSurf& surface,
    const label tI,
    const point& s,
    const point& e,
    point& intersection
)
{
    const pointField& pts = surface.points();
    const labelledTri& tri = surface[tI];

    const triangle<point, point> tria
    (
        pts[tri[0]],
        pts[tri[1]],
        pts[tri[2]]
    );

    return triLineIntersection(tria, s, e, intersection);
}

inline bool boundBoxLineIntersection
(
    const point& s,
    const point& e,
    const boundBox& bb
)
{
    scalar tMax(1.0+SMALL), tMin(-SMALL);

    const vector v = e - s;
    const scalar d = mag(v);

    //- check if the vector has length
    if( d < VSMALL )
    {
        if( bb.contains(s) )
            return true;

        return false;
    }

    const point& pMin = bb.min();
    const point& pMax = bb.max();

    //- check coordinates
    for(label dir=0;dir<3;++dir)
    {
        const scalar vd = v[dir];
        const scalar sd = s[dir];

        if( mag(vd) > (SMALL * d) )
        {
            if( vd >= 0.0 )
            {
                tMin = Foam::max(tMin, (pMin[dir] - sd) / vd);
                tMax = Foam::min(tMax, (pMax[dir] - sd) / vd);
            }
            else
            {
                tMin = Foam::max(tMin, (pMax[dir] - sd) / vd);
                tMax = Foam::min(tMax, (pMin[dir] - sd) / vd);
            }
        }
        else if( (sd < pMin[dir]) || (sd > pMax[dir]) )
        {
            return false;
        }
    }

    if( (tMax - tMin) > -SMALL )
        return true;

    return false;
}

inline bool doFaceAndTriangleIntersect
(
    const triSurf& surface,
    const label triI,
    const face& f,
    const pointField& facePoints
)
{
    const pointField& triPoints = surface.points();

    const point centre = f.centre(facePoints);
    point intersection;

    //- check if any triangle edge intersects the face
    const labelledTri& tri = surface[triI];

    forAll(tri, eI)
    {
        const point& s = triPoints[tri[eI]];
        const point& e = triPoints[tri[(eI+1)%3]];

        forAll(f, pI)
        {
            const triangle<point, point> tria
            (
                facePoints[f[pI]],
                facePoints[f.nextLabel(pI)],
                centre
            );

            const bool currIntersection =
                help::triLineIntersection
                (
                    tria,
                    s,
                    e,
                    intersection
                );

            if( currIntersection )
                return true;
        }
    }

    //- check if any face edges intersect the triangle
    forAll(f, pI)
    {
        const point& s = facePoints[f[pI]];
        const point& e = facePoints[f.nextLabel(pI)];

        const bool intersected =
            help::triLineIntersection
            (
                surface,
                triI,
                s,
                e,
                intersection
            );

        if( intersected )
            return true;
    }

    return false;
}

inline bool pointInsideFace
(
    const point& p,
    const face& f,
    const vector& n,
    const pointField& fp
)
{
    const edgeList fe = f.edges();
    forAll(f, pI)
    {
        if( mag(p - fp[f[pI]]) < SMALL )
            return true;

        vector pv = p - fp[f[pI]];
        pv /= mag(pv);

        vector lv = n ^ fe[pI].vec(fp);
        lv /= mag(lv);

        const scalar d = pv & lv;
        if( d < -SMALL )
        {
            return false;
        }
    }

    return true;
}

inline bool pointInsideFace
(
    const point& p,
    const face& f,
    const pointField& fp
)
{
    vector n = f.normal(fp);
    n /= mag(n);

    const edgeList fe = f.edges();
    forAll(f, pI)
    {
        if( mag(p - fp[f[pI]]) < SMALL )
        return true;

        vector pv = p - fp[f[pI]];
        pv /= mag(pv);

        vector lv = n ^ fe[pI].vec(fp);
        lv /= mag(lv);

        const scalar d = pv & lv;
        if( d < -SMALL )
        {
            return false;
        }
    }

    return true;
}

inline point nearestPointOnTheEdge
(
    const point& edgePoint0,
    const point& edgePoint1,
    const point& p
)
{
    const vector e = edgePoint1 - edgePoint0;
    const scalar d = mag(e);
    const vector k = p - edgePoint0;

    if( d < VSMALL )
        return edgePoint0;

    return edgePoint0 + ((e / (d*d)) * (e & k));
}

inline point nearestPointOnTheEdgeExact
(
    const point& edgePoint0,
    const point& edgePoint1,
    const point& p
)
{
    const vector e = edgePoint1 - edgePoint0;
    const scalar d = mag(e);
    const vector k = p - edgePoint0;

    if( d < VSMALL )
        return edgePoint0;

    const scalar t = (e & k) / (d * d + VSMALL);
    if( t > 1.0 )
    {
        return edgePoint1;
    }
    else if( t < 0.0 )
    {
        return edgePoint0;
    }

    return edgePoint0 + (e * t);
}

inline scalar distanceOfPointFromTheEdge
(
    const point& edgePoint0,
    const point& edgePoint1,
    const point& p
)
{
    return mag(nearestPointOnTheEdge(edgePoint0, edgePoint1, p) - p);
}

inline label numberOfFaceGroups
(
    const labelHashSet& containedElements,
    const point& centre,
    const scalar range,
    const triSurf& surface
)
{
    const pointField& points = surface.points();
    const edgeLongList& edges = surface.edges();
    const VRWGraph& faceEdges = surface.facetEdges();
    const VRWGraph& edgeFaces = surface.edgeFacets();

    labelHashSet triaInRange(containedElements.size());
    const scalar rangeSq = range * range;
    forAllConstIter(labelHashSet, containedElements, it)
    {
        const point p = nearestPointOnTheTriangle(it.key(), surface, centre);
        if( magSqr(p - centre) < rangeSq )
            triaInRange.insert(it.key());
    }

    Map<label> elGroup(triaInRange.size());
    Map<label> testEdge(triaInRange.size());

    label nGroups(0);

    DynList<label> front;

    forAllConstIter(labelHashSet, triaInRange, it)
        if( !elGroup.found(it.key()) )
        {
            front.clear();
            front.append(it.key());
            elGroup.insert(it.key(), nGroups);

            while( front.size() )
            {
                const label fLabel = front.removeLastElement();

                forAllRow(faceEdges, fLabel, feI)
                {
                    const label edgeI = faceEdges(fLabel, feI);

                    //- check if the edge intersects the bounding box
                    if( testEdge.found(edgeI) )
                    {
                        if( !testEdge[edgeI] )
                            continue;
                    }
                    else
                    {
                        const point& s = points[edges[edgeI][0]];
                        const point& e = points[edges[edgeI][1]];
                        const point np = nearestPointOnTheEdge(s, e, centre);
                        if( magSqr(np - centre) < rangeSq )
                        {
                            testEdge.insert(edgeI, 1);
                        }
                        else
                        {
                            testEdge.insert(edgeI, 0);
                            continue;
                        }
                    }

                    forAllRow(edgeFaces, edgeI, efI)
                    {
                        const label nei = edgeFaces(edgeI, efI);
                        if
                        (
                            triaInRange.found(nei) &&
                            !elGroup.found(nei)
                        )
                        {
                            elGroup.insert(nei, nGroups);
                            front.append(nei);
                        }
                    }
                }
            }

            ++nGroups;
        }

    return nGroups;
}

inline label numberOfEdgeGroups
(
    const labelHashSet& containedEdges,
    const point& centre,
    const scalar range,
    const triSurf& surface
)
{
    const pointField& points = surface.points();
    const edgeLongList& edges = surface.edges();
    const VRWGraph& pointEdges = surface.pointEdges();

    const scalar rangeSq = range * range;
    labelHashSet edgesInRange(containedEdges.size());
    forAllConstIter(labelHashSet, containedEdges, it)
    {
        const edge& e = edges[it.key()];
        const point& sp = points[e[0]];
        const point& ep = points[e[1]];

        const point p = nearestPointOnTheEdgeExact(sp, ep, centre);
        if( magSqr(p - centre) < rangeSq )
            edgesInRange.insert(it.key());
    }

    Map<label> elGroup(edgesInRange.size());
    Map<label> pointTest(edgesInRange.size());

    label nGroups(0);

    DynList<label> front;

    forAllConstIter(labelHashSet, edgesInRange, it)
        if( !elGroup.found(it.key()) )
        {
            front.clear();
            front.append(it.key());
            elGroup.insert(it.key(), nGroups);

            while( front.size() )
            {
                const label eLabel = front.removeLastElement();
                const edge& e = edges[eLabel];

                for(label i=0;i<2;++i)
                {
                    if( pointTest.found(e[i]) )
                    {
                        if( !pointTest[e[i]] )
                            continue;
                    }
                    else
                    {
                        if( magSqr(points[e[i]] - centre) < rangeSq )
                        {
                            pointTest.insert(e[i], 1);
                        }
                        else
                        {
                            pointTest.insert(e[i], 0);
                            continue;
                        }
                    }

                    forAllRow(pointEdges, e[i], peI)
                    {
                        const label nei = pointEdges(e[i], peI);

                        if( edgesInRange.found(nei) && !elGroup.found(nei) )
                        {
                            elGroup.insert(nei, nGroups);
                            front.append(nei);
                        }
                    }
                }
            }

            ++nGroups;
        }

    return nGroups;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
