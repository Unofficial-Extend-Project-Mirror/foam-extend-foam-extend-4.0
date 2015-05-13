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

#include "triSurf.H"
#include "meshOctreeCubeCoordinates.H"
#include "helperFunctions.H"

//#define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Static data

const label meshOctreeCubeCoordinates::edgeNodes_[12][2] =
    {
        //- edges in x-direction
        {0, 1},
        {2, 3},
        {4, 5},
        {6, 7},
        //- edges in y-direction
        {0, 2},
        {1, 3},
        {4, 6},
        {5, 7},
        //- edges in z-direction
        {0, 4},
        {1, 5},
        {2, 6},
        {3, 7}
    };

const label meshOctreeCubeCoordinates::faceNodes_[6][4] =
    {
        {0, 4, 6, 2},
        {1, 3, 7, 5},
        {0, 1, 5, 4},
        {2, 6, 7, 3},
        {0, 2, 3, 1},
        {4, 5, 7, 6}
    };

const label meshOctreeCubeCoordinates::nodeFaces_[8][3] =
    {
        {0, 2, 4},
        {1, 2, 4},
        {0, 3, 4},
        {1, 3, 4},
        {0, 2, 5},
        {1, 2, 5},
        {0, 3, 5},
        {1, 3, 5}
    };

const label meshOctreeCubeCoordinates::faceEdges_[6][4] =
    {
        {8, 6, 10, 4},
        {5, 11, 7, 9},
        {0, 9, 2, 8},
        {10, 3, 11, 1},
        {4, 1, 5, 0},
        {2, 7, 3, 6}
    };

const label meshOctreeCubeCoordinates::edgeFaces_[12][2] =
    {
        {2, 4},
        {3, 4},
        {2, 5},
        {3, 5},
        {0, 4},
        {1, 4},
        {0, 5},
        {1, 5},
        {0, 2},
        {1, 2},
        {0, 3},
        {1, 3}
    };

const label meshOctreeCubeCoordinates::oppositeFace_[6] = {1, 0, 3, 2, 5, 4};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeCubeCoordinates::vertices
(
    const boundBox& rootBox,
    FixedList<point, 8>& vrt
) const
{
    const vector tol = SMALL * (rootBox.max() - rootBox.min());

    point min_, max_;
    cubeBox(rootBox, min_, max_);
    min_ -= tol;
    max_ += tol;

    vrt[0] = point(min_.x(), min_.y(), min_.z());
    vrt[1] = point(max_.x(), min_.y(), min_.z());
    vrt[2] = point(min_.x(), max_.y(), min_.z());
    vrt[3] = point(max_.x(), max_.y(), min_.z());
    vrt[4] = point(min_.x(), min_.y(), max_.z());
    vrt[5] = point(max_.x(), min_.y(), max_.z());
    vrt[6] = point(min_.x(), max_.y(), max_.z());
    vrt[7] = point(max_.x(), max_.y(), max_.z());
}

void meshOctreeCubeCoordinates::edgeVertices
(
    const boundBox& rootBox,
    FixedList<FixedList<point, 2>, 12>& e
) const
{
    FixedList<point, 8> vertices;
    this->vertices(rootBox, vertices);

    for(label i=0;i<12;++i)
    {
        e[i][0] = vertices[edgeNodes_[i][0]];
        e[i][1] = vertices[edgeNodes_[i][1]];
    }
}

bool meshOctreeCubeCoordinates::intersectsTriangle
(
    const triSurf& surface,
    const boundBox& rootBox,
    const label tI
) const
{
    const pointField& points = surface.points();
    const labelledTri& ltri = surface[tI];

    const vector tol = SMALL * (rootBox.max() - rootBox.min());

    //- calculate the bound box of the octree cube
    boundBox cBox;
    cubeBox(rootBox, cBox.min(), cBox.max());
    cBox.min() -= tol;
    cBox.max() += tol;

    //- calculate the bounding box of the triangle
    boundBox tBox;
    tBox.min() = tBox.max() = points[ltri[0]];

    for(label pI=1;pI<3;++pI)
    {
        tBox.max() = Foam::max(points[ltri[pI]], tBox.max());
        tBox.min() = Foam::min(points[ltri[pI]], tBox.min());
    }

    tBox.min() -= tol;
    tBox.max() += tol;

    return cBox.overlaps(tBox);
}

bool meshOctreeCubeCoordinates::intersectsTriangleExact
(
    const triSurf& surface,
    const boundBox& rootBox,
    const label tI
) const
{
    if( !intersectsTriangle(surface, rootBox, tI) )
        return false;

    const vector tol = SMALL * (rootBox.max() - rootBox.min());

    const pointField& points = surface.points();
    const labelledTri& ltri = surface[tI];

    //- check if any of the vertices is in the cube
    forAll(ltri, pI)
        if( isVertexInside(rootBox, points[ltri[pI]]) )
            return true;

    //- check if any edges of the triangle intersect the cube
    boundBox bb;
    cubeBox(rootBox, bb.min(), bb.max());
    bb.min() -= tol;
    bb.max() += tol;

    for(label eI=0;eI<3;++eI)
    {
        const edge edg(ltri[eI], ltri[(eI+1)%3]);
        const point& s = points[edg.start()];
        const point& e = points[edg.end()];

        if( help::boundBoxLineIntersection(s, e, bb) )
            return true;
    }

    //- check if any cube edges intersects the triangle
    FixedList<FixedList<point, 2>, 12> e;
    this->edgeVertices(rootBox, e);

    point intersection;
    forAll(e, eI)
        if(
            help::triLineIntersection
            (
                surface,
                tI,
                e[eI][0],
                e[eI][1],
                intersection
            )
        )
            return true;

    return false;
}

bool meshOctreeCubeCoordinates::isVertexInside
(
    const boundBox& rootBox,
    const point& p
) const
{
    const vector tol = SMALL * (rootBox.max() - rootBox.min());

    point min, max;
    cubeBox(rootBox, min, max);
    max += tol;
    min -= tol;

    if(
        ((p.x() - max.x()) > 0.0) ||
        ((p.y() - max.y()) > 0.0) ||
        ((p.z() - max.z()) > 0.0) ||
        ((p.x() - min.x()) < 0.0) ||
        ((p.y() - min.y()) < 0.0) ||
        ((p.z() - min.z()) < 0.0)
    )
        return false;

    return true;
}

bool meshOctreeCubeCoordinates::isPositionInside
(
    const meshOctreeCubeCoordinates& cc
) const
{
    # ifdef DEBUGSearch
    Info << "Checking cube " << *this << endl;
    Info << "level " << short(l) << " x: " << px
        << " y: " << py << " z:" << pz << endl;
    # endif

    if( cc.level() >= this->level() )
    {
        const direction diff = cc.level() - this->level();
        meshOctreeCubeCoordinates reducedLevel =
            cc.reduceLevelBy(diff);

        # ifdef DEBUGSearch
        Info << "diff " << label(diff) << endl;
        Info << "Divider " << divider << endl;
        Info << "Coordinates at level are " << reducedLevel << endl;
        # endif

        if( reducedLevel == *this )
            return true;
    }
    else
    {
        FatalErrorIn
        (
            "bool meshOctreeCubeCoordinates::isPositionInside"
            "(const label px,const label py,"
            "const label pz,const direction l)"
        ) << "Cannot find exact position of finer cube" << exit(FatalError);
    }

    return false;
}

bool meshOctreeCubeCoordinates::intersectsLine
(
    const boundBox& rootBox,
    const point& s,
    const point& e
) const
{
    const scalar tol = SMALL * (rootBox.max().x() - rootBox.min().x());

    point min, max;
    cubeBox(rootBox, min, max);

    //- check if the cube contains start point or end point
    min -= vector(tol,tol,tol);
    max += vector(tol,tol,tol);

    //- check for intersections of line with the cube faces
    const vector v(e - s);
    scalar t;
    point i;

    if( mag(v.x()) > tol )
    {
        //- x-min face
        t = (min.x() - s.x()) / v.x();
        i = s + t * v;
        if(
            (t > -tol) && (t < (1.0+tol)) &&
            (i.y() - min.y() > -tol) &&
            (i.y() - max.y() < tol) &&
            (i.z() - min.z() > -tol) &&
            (i.z() - max.z() < tol)
        )
            return true;

        //- x-max face
        t = (max.x() - s.x()) / v.x();
        i = s + t * v;
        if(
            (t > -tol) && (t < (1.0+tol)) &&
            (i.y() - min.y() > -tol) &&
            (i.y() - max.y() < tol) &&
            (i.z() - min.z() > -tol) &&
            (i.z() - max.z() < tol)
        )
            return true;
    }

    if( mag(v.y()) > tol)
    {
        //- y-min face
        t = (min.y() - s.y()) / v.y();
        i = s + t * v;
        if(
            (t > -tol) && (t < (1.0+tol)) &&
            (i.x() - min.x() > -tol) &&
            (i.x() - max.x() < tol) &&
            (i.z() - min.z() > -tol) &&
            (i.z() - max.z() < tol)
        )
            return true;

        //- y-max face
        t = (max.y() - s.y()) / v.y();
        i = s + t * v;
        if(
            (t > -tol) && (t < (1.0+tol)) &&
            (i.x() - min.x() > -tol) &&
            (i.x() - max.x() < tol) &&
            (i.z() - min.z() > -tol) &&
            (i.z() - max.z() < tol)
        )
            return true;
    }

    if( mag(v.z()) > tol )
    {
        //- z-min face
        t = (min.z() - s.z()) / v.z();
        i = s + t * v;
        if(
            (t > -tol) && (t < (1.0+tol)) &&
            (i.x() - min.x() > -tol) &&
            (i.x() - max.x() < tol) &&
            (i.y() - min.y() > -tol) &&
            (i.y() - max.y() < tol)
        )
            return true;

        //- z-min face
        t = (max.z() - s.z()) / v.z();
        i = s + t * v;
        if(
            (t > -tol) && (t < (1.0+tol)) &&
            (i.x() - min.x() > -tol) &&
            (i.x() - max.x() < tol) &&
            (i.y() - min.y() > -tol) &&
            (i.y() - max.y() < tol)
        )
            return true;
    }

    if( isVertexInside(rootBox, s) )
        return true;

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
