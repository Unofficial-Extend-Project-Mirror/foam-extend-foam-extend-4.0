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
#include "Pair.H"
#include "Map.H"

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
    const triangle<point, point>& tri,
    const point& p
)
{
    const vector v0 = tri.b() - tri.a();
    const vector v1 = tri.c() - tri.a();
    const vector v2 = p - tri.a();

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

        FixedList<Pair<point>, 3> edges;
        edges[0] = Pair<point>(tri.a(), tri.b());
        edges[1] = Pair<point>(tri.b(), tri.c());
        edges[2] = Pair<point>(tri.c(), tri.a());

        forAll(edges, eI)
        {
            const Pair<point>& e = edges[eI];
            const point np =
                nearestPointOnTheEdge
                (
                    e.first(),
                    e.second(),
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

    const point pProj = tri.a() + u * v0 + v * v1;

    //- Check if point is in triangle
    if( (u >= -SMALL) && (v >= -SMALL) && ((u + v) <= (1.0+SMALL)) )
    {
        return pProj;
    }
    else
    {
        if( u < -SMALL )
        {
            const scalar ed = ((pProj - tri.a()) & v1) / (magSqr(v1) + VSMALL);

            if( ed > 1.0 )
            {
                return tri.c();
            }
            else if( ed < 0.0 )
            {
                return tri.a();
            }
            else
            {
                return (tri.a() + v1 * ed);
            }
        }
        else if( v < -SMALL )
        {
            const scalar ed = ((pProj - tri.a()) & v0) / (magSqr(v0) + VSMALL);

            if( ed > 1.0 )
            {
                return tri.b();
            }
            else if( ed < 0.0 )
            {
                return tri.a();
            }
            else
            {
                return (tri.a() + v0 * ed);
            }
        }
        else
        {
            const vector ev = tri.b() - tri.c();
            const scalar ed = ((pProj - tri.c()) & ev) / (magSqr(ev) + VSMALL);

            if( ed > 1.0 )
            {
                return tri.b();
            }
            else if( ed < 0.0 )
            {
                return tri.c();
            }
            else
            {
                return (tri.c() + ev * ed);
            }
        }
    }

    return pProj;
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

    const triangle<point, point> tri
    (
        points[ltri[0]],
        points[ltri[1]],
        points[ltri[2]]
    );

    return nearestPointOnTheTriangle(tri, p);
}

inline bool findMinimizerPoint
(
    const DynList<point>& origins,
    const DynList<vector>& normals,
    point& pMin
)
{
    if( origins.size() != normals.size() )
        FatalErrorIn
        (
            "inline bool findMinimizerPoint"
            "(const DynList<point>&, const DynList<vector>&, point&)"
        ) << "Size of normals and origins do not match" << abort(FatalError);

    tensor mat(tensor::zero);
    vector source(vector::zero);

    forAll(origins, i)
    {
        //- use the normalized vector
        vector n = normals[i];
        n /= (mag(n) + VSMALL);

        const tensor t = n * n;

        mat += t;

        source += (origins[i] & n) * n;
    }

    const scalar determinant = Foam::mag(Foam::det(mat));
    if( determinant < SMALL )
        return false;

    pMin = (inv(mat, determinant) & source);

    return true;
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

inline bool lineFaceIntersection
(
    const point& sp,
    const point& ep,
    const face& f,
    const pointField& fp,
    point& intersection
)
{
    const point c = f.centre(fp);

    forAll(f, pI)
    {
        const triangle<point, point> tri
        (
            fp[f[pI]],
            fp[f.nextLabel(pI)],
            c
        );

        if( triLineIntersection(tri, sp, ep, intersection) )
            return true;
    }

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

inline bool doEdgesOverlap
(
    const point& e0p0,
    const point& e0p1,
    const point& e1p0,
    const point& e1p1,
    FixedList<point, 2>& overlappingPart,
    const scalar distTol,
    const scalar cosTol
)
{
    if( distTol < 0.0 )
    {
        WarningIn
        (
            "inline bool doEdgesOverlap(const point, const point&,"
            "const point&, const point&, const scalar, const scalar,"
            "FixedList<point, 2>& overlappingPart)"
        ) << "Distance is not specified" << endl;

        return false;
    }

    vector e0 = e0p1 - e0p0;
    const scalar de0 = mag(e0);
    e0 /= (de0 + VSMALL);

    vector e1 = e1p1 - e1p0;
    const scalar de1 = mag(e1);
    e1 /= (de1 + VSMALL);

    //- check the angle deviation between the two vectors
    if( mag(e0 & e1) < cosTol )
        return false;

    scalar t00 = (e1p0 - e0p0) & e0;
    scalar t01 = (e1p1 - e0p0) & e0;

    scalar t10 = (e0p0 - e1p0) & e1;
    scalar t11 = (e0p1 - e1p0) & e1;

    //- check if points are colinear within the tolerance
    if
    (
        (magSqr(e0p0 + t00*e0 - e1p0) <= distTol*distTol) ||
        (magSqr(e0p0 + t01*e0 - e1p1) <= distTol*distTol) ||
        (magSqr(e1p0 + t10*e1 - e0p0) <= distTol*distTol) ||
        (magSqr(e1p0 + t11*e1 - e0p1) <= distTol*distTol)
    )
    {
        vector vec = e0 + (((e0 & e1) > 0.0) ? e1 : -e1);
        vec /= (mag(vec) + VSMALL);
        const point origin = 0.25 * (e0p0 + e0p1 + e1p0 + e1p1);

        //- calculate parameters on the regression line
        t00 = (e0p0 - origin) & vec;
        t01 = (e0p1 - origin) & vec;
        t10 = (e1p0 - origin) & vec;
        t11 = (e1p1 - origin) & vec;

        //- check interval overlapping over the line
        const scalar t0Min = Foam::min(t00, t01);
        const scalar t0Max = Foam::max(t00, t01);
        const scalar t1Min = Foam::min(t10, t11);
        const scalar t1Max = Foam::max(t10, t11);

        if( t1Min < t0Max )
        {
            overlappingPart[0] = origin + t1Min * vec;
            overlappingPart[1] = origin + t0Max * vec;

            return true;
        }
        else if( t0Min < t1Max )
        {
            overlappingPart[0] = origin + t0Min * vec;
            overlappingPart[1] = origin + t1Max * vec;

            return true;
        }
    }

    return false;
}

inline bool doTrianglesOverlap
(
    const triangle<point, point>& tri0,
    const triangle<point, point>& tri1,
    DynList<point>& overlappingPolygon,
    const scalar distTol,
    const scalar cosTol
)
{
    if( distTol < 0.0 )
    {
        WarningIn
        (
            "inline bool doTrianglesOverlap(const triangle<point, point>&,"
            " const triangle<point, point>&, DynList<point>&,"
            " const scalar, const scalar)"
        ) << "Distance is not specified" << endl;

        return false;
    }

    vector n0 = tri0.normal();
    const scalar dn0 = mag(n0);
    n0 /= (dn0 + VSMALL);

    if( dn0 < VSMALL )
        return false;

    vector n1 = tri1.normal();
    const scalar dn1 = mag(n1);
    n1 /= (dn1 + VSMALL);

    if( dn1 < VSMALL )
        return false;

    //- check the angle deviation between the two vectors
    if( (mag(n0 & n1) < cosTol) && (dn0 >= VSMALL) && (dn1 >= VSMALL) )
        return false;

    //- check if the two nearest points are within tolerance
    if
    (
        (mag((tri1.a() - tri0.a()) & n0) < distTol) ||
        (mag((tri1.b() - tri0.a()) & n0) < distTol) ||
        (mag((tri1.c() - tri0.a()) & n0) < distTol) ||
        (mag((tri0.a() - tri1.a()) & n1) < distTol) ||
        (mag((tri0.b() - tri1.a()) & n1) < distTol) ||
        (mag((tri0.c() - tri1.a()) & n1) < distTol)
    )
    {
        vector vec = n0 + ((n0 & n1) >= 0.0 ? n1 : -n1);
        vec /= (mag(vec) + VSMALL);
        const point origin = 0.5 * (tri0.centre() + tri1.centre());

        overlappingPolygon.clear();

        //- calculate max distance point from the origin
        point bestPoint = tri0.a();
        scalar distSq = magSqr(tri0.a() - origin);

        scalar dSq = magSqr(tri0.b() - origin);
        if( dSq > distSq )
        {
            distSq = dSq;
            bestPoint = tri0.b();
        }

        dSq = magSqr(tri0.c() - origin);
        if( dSq > distSq )
        {
            distSq = dSq;
            bestPoint = tri0.c();
        }

        dSq = magSqr(tri1.a() - origin);
        if( dSq > distSq )
        {
            distSq = dSq;
            bestPoint = tri1.a();
        }

        dSq = magSqr(tri1.b() - origin);
        if( dSq > distSq )
        {
            distSq = dSq;
            bestPoint = tri1.b();
        }

        dSq = magSqr(tri1.c() - origin);
        if( dSq > distSq )
        {
            distSq = dSq;
            bestPoint = tri1.c();
        }

        if( distSq < VSMALL )
            return false;

        //- transform into planar coordinates
        vector x = (bestPoint - origin) - ((bestPoint - origin) & vec) * vec;
        x /= (mag(x) + VSMALL);
        vector y = vec ^ x;

        DynList<point2D, 6> poly2D(3);
        poly2D[0] = point2D((tri0.a() - origin) & x, (tri0.a() - origin) & y);
        poly2D[1] = point2D((tri0.b() - origin) & x, (tri0.b() - origin) & y);
        poly2D[2] = point2D((tri0.c() - origin) & x, (tri0.c() - origin) & y);

        FixedList<point2D, 3> t1Proj;
        t1Proj[0] = point2D((tri1.a() - origin) & x, (tri1.a() - origin) & y);
        t1Proj[1] = point2D((tri1.b() - origin) & x, (tri1.b() - origin) & y);
        t1Proj[2] = point2D((tri1.c() - origin) & x, (tri1.c() - origin) & y);

        forAll(t1Proj, eI)
        {
            const vector2D vec = t1Proj[(eI+1)%3] - t1Proj[eI];

            DynList<scalar, 6> distance(poly2D.size());
            forAll(poly2D, pI)
            {
                const vector2D pVec = poly2D[pI] - t1Proj[eI];
                distance[pI] = vec.y() * pVec.x() - vec.x() * pVec.y();
            }

            DynList<point2D, 6> newPoly2D;

            forAll(distance, pI)
            {
                if( distance[pI] >= 0.0 )
                {
                    newPoly2D.append(poly2D[pI]);

                    if( distance.fcElement(pI) < 0.0 )
                    {
                        //- this is very sensitive to floaing point tolerances
                        const point2D newP =
                            (
                                mag(distance[pI]) * poly2D.fcElement(pI) +
                                mag(distance.fcElement(pI)) * poly2D[pI]
                            ) /
                            (
                                mag(distance[pI]) +
                                mag(distance.fcElement(pI)) +
                                VSMALL
                            );

                        newPoly2D.append(newP);
                    }
                }
                else if( distance.fcElement(pI) >= 0.0 )
                {
                    //- this is very sensitive to floaing point tolerances
                    const point2D newP =
                        (
                            mag(distance[pI]) * poly2D.fcElement(pI) +
                            mag(distance.fcElement(pI)) * poly2D[pI]
                        ) /
                        (
                            mag(distance[pI]) +
                            mag(distance.fcElement(pI)) +
                            VSMALL
                        );

                    newPoly2D.append(newP);
                }
            }

            poly2D = newPoly2D;
        }

        //- check if the overlapping polygon exists
        if( poly2D.size() == 0 )
        {
            overlappingPolygon.clear();
            return false;
        }

        //- fill the overlapping polygon
        overlappingPolygon.setSize(poly2D.size());
        forAll(poly2D, pI)
        {
            const point2D& pp = poly2D[pI];
            overlappingPolygon[pI] = origin + x * pp.x() + y * pp.y();
        }

        return true;
    }

    overlappingPolygon.clear();
    return false;
}

inline bool doTrianglesIntersect
(
    const triangle<point, point> &tri0,
    const triangle<point, point> &tri1,
    const scalar distTol
)
{
    if( distTol < 0.0 )
    {
        WarningIn
        (
            "inline bool doTrianglesIntersect(const triangle<point, point>&,"
            " const triangle<point, point>&,"
            " const scalar, const scalar)"
        ) << "Distance is not specified" << endl;

        return false;
    }

    //- find distances of points from the second triangle
    point np = nearestPointOnTheTriangle(tri1, tri0.a());
    if( magSqr(np - tri0.a()) < distTol*distTol )
        return true;
    np = nearestPointOnTheTriangle(tri1, tri0.b());
    if( magSqr(np - tri0.b()) < distTol*distTol )
        return true;
    np = nearestPointOnTheTriangle(tri1, tri0.c());
    if( magSqr(np - tri0.c()) < distTol*distTol )
        return true;

    //- find distances of points from the first triangle
    np = nearestPointOnTheTriangle(tri0, tri1.a());
    if( magSqr(np - tri1.a()) < distTol*distTol )
        return true;
    np = nearestPointOnTheTriangle(tri0, tri1.b());
    if( magSqr(np - tri1.b()) < distTol*distTol )
        return true;
    np = nearestPointOnTheTriangle(tri0, tri1.c());
    if( magSqr(np - tri1.c()) < distTol*distTol )
        return true;

    //- find edge intersections
    point intersection;
    if( triLineIntersection(tri0, tri1.a(), tri1.b(), intersection) )
        return true;
    if( triLineIntersection(tri0, tri1.b(), tri1.c(), intersection) )
        return true;
    if( triLineIntersection(tri0, tri1.c(), tri1.a(), intersection) )
        return true;

    if( triLineIntersection(tri1, tri0.a(), tri0.b(), intersection) )
        return true;
    if( triLineIntersection(tri1, tri0.b(), tri0.c(), intersection) )
        return true;
    if( triLineIntersection(tri1, tri0.c(), tri0.a(), intersection) )
        return true;

    return false;
}

inline bool doTrianglesIntersect
(
    const triangle<point, point>& tri0,
    const triangle<point, point>& tri1,
    DynList<point> &intersectionPoints,
    const scalar distTol
)
{
    if( distTol < 0.0 )
    {
        WarningIn
        (
            "inline bool doTrianglesIntersect(const triangle<point, point>&,"
            " const triangle<point, point>&, DynList<point>&,"
            " const scalar, const scalar)"
        ) << "Distance is not specified" << endl;

        return false;
    }

    vector n0 = tri0.normal();
    const scalar dn0 = mag(n0);
    n0 /= (dn0 + VSMALL);

    vector n1 = tri1.normal();
    const scalar dn1 = mag(n1);
    n1 /= (dn1 + VSMALL);

    //- distance of the points of the first triangle from the plane
    //- of the second triangle
    FixedList<scalar, 3> distancesTri0;
    distancesTri0[0] = (tri0.a() - tri1.a()) & n1;
    distancesTri0[1] = (tri0.b() - tri1.a()) & n1;
    distancesTri0[2] = (tri0.c() - tri1.a()) & n1;

    forAll(distancesTri0, i)
        if( mag(distancesTri0[i]) < distTol )
            distancesTri0[i] = 0.0;

    bool hasPositive(false), hasNegative(false), hasZero(false);
    DynList<point, 2> intersectionLine0;
    forAll(distancesTri0, pI)
    {
        if( distancesTri0[pI] >= distTol )
        {
            hasPositive = true;
        }
        else if( distancesTri0[pI] <= -distTol )
        {
            hasNegative = true;
        }
        else
        {
            hasZero = true;
        }
    }

    //- find points on the intersection line
    if
    (
        (hasPositive && !hasNegative && !hasZero) ||
        (hasNegative && !hasPositive && !hasZero)
    )
    {
        //- there can be no intersection
        intersectionPoints.clear();
        return false;
    }
    else if
    (
        (hasPositive && (hasNegative || hasZero)) ||
        (hasNegative && (hasPositive || hasZero))
    )
    {
        //- find points on the intersection line
        if( mag(distancesTri0[0]) < distTol)
        {
            intersectionLine0.append(tri0.a());
        }
        else if( distancesTri0[0] * distancesTri0[1] < 0.0 )
        {
            intersectionLine0.append
            (
                (tri0.a() * distancesTri0[1] - tri0.b() * distancesTri0[0]) /
                (distancesTri0[0] - distancesTri0[1])
            );
        }

        if( mag(distancesTri0[1]) < distTol )
        {
            intersectionLine0.append(tri0.b());
        }
        else if( distancesTri0[1] * distancesTri0[2] < 0.0 )
        {
            intersectionLine0.append
            (
                (tri0.b() * distancesTri0[2] - tri0.c() * distancesTri0[1]) /
                (distancesTri0[1] - distancesTri0[2])
            );
        }

        if( mag(distancesTri0[2]) < distTol )
        {
            intersectionLine0.append(tri0.c());
        }
        else if( distancesTri0[2] * distancesTri0[0] < 0.0 )
        {
            intersectionLine0.append
            (
                (tri0.c() * distancesTri0[0] - tri0.a() * distancesTri0[2]) /
                (distancesTri0[2] - distancesTri0[0])
            );
        }
    }
    else
    {
        //- these triangles are in the same plane
        //- check if they overlap
        return doTrianglesOverlap(tri0, tri1, intersectionPoints, distTol);
    }

    //- distance of the points of the second triangle from the plane
    //- of the first triangle
    FixedList<scalar, 3> distancesTri1;
    distancesTri1[0] = (tri1.a() - tri0.a()) & n0;
    distancesTri1[1] = (tri1.b() - tri0.a()) & n0;
    distancesTri1[2] = (tri1.c() - tri0.a()) & n0;

    forAll(distancesTri1, i)
        if( mag(distancesTri1[i]) < distTol )
            distancesTri1[i] = 0.0;

    hasPositive = false;
    hasNegative = false;
    hasZero = false;

    DynList<point, 2> intersectionLine1;
    forAll(distancesTri1, pI)
    {
        if( distancesTri1[pI] >= distTol )
        {
            hasPositive = true;
        }
        else if( distancesTri1[pI] <= -distTol )
        {
            hasNegative = true;
        }
        else
        {
            hasZero = true;
        }
    }

    if
    (
        (hasPositive && !hasNegative && !hasZero) ||
        (hasNegative && !hasPositive && !hasZero)
    )
    {
        //- there can be no intersection
        intersectionPoints.clear();
        return false;
    }
    else if
    (
        (hasPositive && (hasNegative || hasZero)) ||
        (hasNegative && (hasPositive || hasZero))
    )
    {
        //- find points on the intersection line
        if( mag(distancesTri1[0]) < distTol)
        {
            intersectionLine1.append(tri1.a());
        }
        else if( distancesTri1[0] * distancesTri1[1] < 0.0 )
        {
            intersectionLine1.append
            (
                (tri1.a() * distancesTri1[1] - tri1.b() * distancesTri1[0]) /
                (distancesTri1[0] - distancesTri1[1])
            );
        }

        if( mag(distancesTri1[1]) < distTol)
        {
            intersectionLine1.append(tri1.b());
        }
        else if( distancesTri1[1] * distancesTri1[2] < 0.0 )
        {
            intersectionLine1.append
            (
                (tri1.b() * distancesTri1[2] - tri1.c() * distancesTri1[1]) /
                (distancesTri1[1] - distancesTri1[2])
            );
        }

        if( mag(distancesTri1[2]) < distTol)
        {
            intersectionLine1.append(tri1.c());
        }
        else if( distancesTri1[2] * distancesTri1[0] < 0.0 )
            intersectionLine1.append
            (
                (tri1.c() * distancesTri1[0] - tri1.a() * distancesTri1[2]) /
                (distancesTri1[2] - distancesTri1[0])
            );
    }
    else
    {
        //- these triangles are in the same plane
        //- check if they overlap
        return doTrianglesOverlap(tri0, tri1, intersectionPoints, distTol);
    }

    //- calculate the direction of the intersection line
    vector vec = n0 ^ n1;

    //- n0 and n1 are nearly coplanar, try overlap intersection
    if( magSqr(vec) < SMALL )
    {
        //- these triangles are in the same plane
        return doTrianglesOverlap(tri0, tri1, intersectionPoints, distTol);
    }

    //- calculate intervals on the intersection line
    scalar t0Min(VGREAT), t0Max(-VGREAT);
    point p0Min(vector::zero), p0Max(vector::zero);
    forAll(intersectionLine0, i)
    {
        const scalar t = intersectionLine0[i] & vec;

        if( t < t0Min )
        {
            t0Min = t;
            p0Min = intersectionLine0[i];
        }
        if( t > t0Max )
        {
            t0Max = t;
            p0Max = intersectionLine0[i];
        }
    }

    scalar t1Min(VGREAT), t1Max(-VGREAT);
    point p1Min(vector::zero), p1Max(vector::zero);
    forAll(intersectionLine1, i)
    {
        const scalar t = intersectionLine1[i] & vec;

        if( t < t1Min )
        {
            t1Min = t;
            p1Min = intersectionLine1[i];
        }
        if( t > t1Max )
        {
            t1Max = t;
            p1Max = intersectionLine1[i];
        }
    }

    if( (t1Min <= t0Max) && (t1Max >= t0Min) )
    {
        intersectionPoints.setSize(2);
        intersectionPoints[0] = p1Min;
        intersectionPoints[1] = p0Max;

        return true;
    }
    else if( (t0Min <= t1Max) && (t0Max >= t1Min) )
    {
        intersectionPoints.setSize(2);
        intersectionPoints[0] = p0Min;
        intersectionPoints[1] = p1Max;

        return true;
    }

    intersectionPoints.setSize(0);
    return false;
}

inline bool pointInsideFace
(
    const point& p,
    const face& f,
    const vector& n,
    const pointField& fp,
    const scalar distTol
)
{
    const edgeList fe = f.edges();
    forAll(f, pI)
    {
        if( mag(p - fp[f[pI]]) < distTol )
            return true;

        vector pv = p - fp[f[pI]];
        pv /= mag(pv);

        vector lv = n ^ fe[pI].vec(fp);
        lv /= mag(lv);

        const scalar d = pv & lv;
        if( d < -distTol )
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
    const pointField& fp,
    const scalar distTol
)
{
    const point c = f.centre(fp);
    const scalar tolSqr = sqr(distTol);

    forAll(f, eI)
    {
        const edge fe = f.faceEdge(eI);

        triangle<point, point> tri
        (
            fp[fe[0]],
            fp[fe[1]],
            c
        );

        const point np = nearestPointOnTheTriangle(tri, p);

        if( magSqr(np - p) <= tolSqr )
            return true;
    }

    return false;
}

inline bool isFaceConvexAndOk
(
    const face& f,
    const pointField& fp,
    DynList<bool>& OkPoints
)
{
    bool valid(true);

    vector normal = f.normal(fp);
    const scalar magN = mag(normal);

    //- face has zero area. All points are inverted
    if( magN < VSMALL )
    {
        OkPoints.setSize(f.size());
        OkPoints = false;

        return false;
    }

    normal /= (magN + VSMALL);

    //- project points into the plane formed by the face centre
    //- and the face normal
    plane pl(f.centre(fp), normal);

    DynList<point> projFace(f.size());
    forAll(f, pI)
        projFace[pI] = pl.nearestPoint(fp[f[pI]]);

    //- calculate normalised edge vectors
    DynList<vector> edgeVecs(f.size());
    forAll(f, eI)
    {
        vector& v = edgeVecs[eI];

        const edge e = f.faceEdge(eI);
        v = e.vec(fp);
        v /= (mag(v) + VSMALL);
    }

    //- check if the face is convex
    OkPoints.setSize(f.size());
    forAll(f, pI)
    {
        const vector& pv = edgeVecs[f.rcIndex(pI)];
        const vector& nv = edgeVecs[pI];

        if( ((pv ^ nv) & normal) >= -0.05 )
        {
            OkPoints[pI] = true;
        }
        else
        {
            OkPoints[pI] = false;
            valid = false;
        }
    }

    return valid;
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

    if( d < ROOTVSMALL )
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
    const scalar dSq = magSqr(e);
    const vector k = p - edgePoint0;

    if( dSq < VSMALL )
        return edgePoint0;

    const scalar t = (e & k) / (dSq + VSMALL);
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
