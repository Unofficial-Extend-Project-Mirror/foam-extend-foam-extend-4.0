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

InNamespace
    Foam::help

Description
    Geometry queries useful for mesh generation

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef helperFunctionsGeometryQueries_H
#define helperFunctionsGeometryQueries_H

#include "DynList.H"
#include "plane.H"
#include "face.H"

#include "triSurf.H"
#include "triangle.H"
#include "tetrahedron.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class boolList;
class labelHashSet;
class boundBox;

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Namespace help functions Declaration
\*---------------------------------------------------------------------------*/

namespace help
{
    //- check if a list has nan entries
    template<class ListType>
    bool isnan(const ListType&);

    //- check if a list has inf entries
    template<class ListType>
    bool isinf(const ListType&);

    //- check if the faces share a convex edge
    template<class Face1, class Face2>
    inline bool isSharedEdgeConvex
    (
        const pointField& points,
        const Face1& f1,
        const Face2& f2
    );

    //- angle between the two faces in radians
    template<class Face1, class Face2>
    inline scalar angleBetweenFaces
    (
        const pointField& points,
        const Face1& f1,
        const Face2& f2
    );

    // merge faces which are in the same patch
    // This is need to merge the triagles which are generated for cells
    // cut by more than on boundary regions
    inline faceList mergePatchFaces
    (
        const List< DynList<label> >& pfcs,
        const pointField& polyPoints
    );

    //- check if the point p belongs to the edge e
    inline bool vertexOnLine
    (
        const point& p,
        const edge& e,
        const pointField& ep
    );

    //- check if the point p belongs to the plane
    inline bool vertexInPlane(const point& p, const plane& pl);

    //- find the vertex on the line of the edge nearest to the point p
    inline point nearestPointOnTheEdge
    (
        const point& edgePoint0,
        const point& edgePoint1,
        const point& p
    );

    //- find the vertex on the edge nearest to the point p
    inline point nearestPointOnTheEdgeExact
    (
        const point& edgePoint0,
        const point& edgePoint1,
        const point& p
    );

    //- find and return the distance between the edge and the point p
    inline scalar distanceOfPointFromTheEdge
    (
        const point& edgePoint0,
        const point& edgePoint1,
        const point& p
    );

    //- find the nearest points on the edge and the line
    inline bool nearestEdgePointToTheLine
    (
        const point& edgePoint0,
        const point& edgePoint1,
        const point& lp0,
        const point& lp1,
        point& nearestOnEdge,
        point& nearestOnLine
    );

    //- check if the edge intersects the plane
    inline bool planeIntersectsEdge
    (
        const point& start,
        const point& end,
        const plane& pl,
        point& intersection
    );

    //- check if a vertex lies inside the tetrahedron
    inline bool pointInTetrahedron
    (
        const point& p,
        const tetrahedron<point, point>& tet
    );

    //- check if a line intersects the triangle, and return the intersection
    inline bool triLineIntersection
    (
        const triangle<point, point>& tria,
        const point& lineStart,
        const point& lineEnd,
        point& intersection
    );

    //- check if a line intersects the triangle and return the intersection
    inline bool triLineIntersection
    (
        const triSurf&,
        const label,
        const point&,
        const point&,
        point&
    );

    //- check if the line intersects the bounding box
    inline bool boundBoxLineIntersection
    (
        const point&,
        const point&,
        const boundBox&
    );

    //- check if the surface triangle and the face intersect
    inline bool doFaceAndTriangleIntersect
    (
        const triSurf& surface,
        const label triI,
        const face& f,
        const pointField& facePoints
    );

    //- find the nearest vertex on the surface triangle to the given point
    inline point nearestPointOnTheTriangle
    (
        const label,
        const triSurf&,
        const point&
    );

    //- check if the point is inside or outside the face
    inline bool pointInsideFace
    (
        const point& p,
        const face& f,
        const vector& n,
        const pointField& fp
    );

    //- check if the point is inside or outside the face
    inline bool pointInsideFace
    (
        const point& p,
        const face& f,
        const pointField& fp
    );

    //- check if the vertex is on the positive side of the face plane
    inline bool isVertexVisible(const point& p, const plane& pl);

    //- find number of face groups within a given range
    inline label numberOfFaceGroups
    (
        const labelHashSet& containedElements,
        const point& centre,
        const scalar range,
        const triSurf& surface
    );

    //- find the number of edge groups within the given range
    inline label numberOfEdgeGroups
    (
        const labelHashSet& containedEdges,
        const point& centre,
        const scalar range,
        const triSurf& surface
    );

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "helperFunctionsGeometryQueriesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
