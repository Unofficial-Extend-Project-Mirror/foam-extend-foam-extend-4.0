/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "triangleFuncs.H"
#include "pointField.H"
#include "treeBoundBox.H"
#include "SortableList.H"
#include "boolList.H"
#include "scalar.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::triangleFuncs::setIntersection
(
    const point& oppositeSidePt,
    const scalar oppositeSign,

    const point& thisSidePt,
    const scalar thisSign,

    const scalar tol,

    point& pt
)
{
    scalar denom = oppositeSign - thisSign;

    if (mag(denom) < tol)
    {
        // If almost does not cut choose one which certainly cuts.
        pt = oppositeSidePt;
    }
    else
    {
        pt = oppositeSidePt + oppositeSign/denom*(thisSidePt - oppositeSidePt);
    }
}


void Foam::triangleFuncs::selectPt
(
    const bool select0,
    const point& p0,
    const point& p1,
    point& min
)
{
    if (select0)
    {
        min = p0;
    }
    else
    {
        min = p1;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Intersect triangle with parallel edges aligned with axis i0.
// Returns true (and intersection in pInter) if any of them intersects triangle.
bool Foam::triangleFuncs::intersectAxesBundle
(
    const point& V0,
    const point& V10,
    const point& V20,
    const label i0,
    const pointField& origin,
    const scalar maxLength,
    point& pInter
)
{
    // Based on Graphics Gems - Fast Ray Triangle intersection.
    // Since direction is coordinate axis there is no need to do projection,
    // we can directly check u,v components for inclusion in triangle.

    scalar localScale = max(max(magSqr(V10), magSqr(V20)), 1.0);

    // Get other components
    label i1 = (i0 + 1) % 3;
    label i2 = (i1 + 1) % 3;

    scalar u1 = V10[i1];
    scalar v1 = V10[i2];

    scalar u2 = V20[i1];
    scalar v2 = V20[i2];

    scalar det = v2*u1 - u2*v1;

    // Fix for  V0:(-31.71428 0 -15.10714)
    //          V10:(-1.285715 8.99165e-16 -1.142858)
    //          V20:(0 0 -1.678573)
    //          i0:0
    if (Foam::mag(det)/localScale < SMALL)
    {
        // Triangle parallel to dir
        return false;
    }

    forAll(origin, originI)
    {
        const point& P = origin[originI];

        scalar u0 = P[i1] - V0[i1];
        scalar v0 = P[i2] - V0[i2];

        scalar alpha = 0;
        scalar beta = 0;
        bool inter = false;

        if (Foam::mag(u1) < ROOTVSMALL)
        {
            beta = u0/u2;
            if ((beta >= 0) && (beta <= 1))
            {
                alpha = (v0 - beta*v2)/v1;
                inter = ((alpha >= 0) && ((alpha + beta) <= 1));
            }
        }
        else
        {
            beta = (v0*u1 - u0*v1)/det;
            if ((beta >= 0) && (beta <= 1))
            {
                alpha = (u0 - beta*u2)/u1;
                inter = ((alpha >= 0) && ((alpha + beta) <= 1));
            }
        }

        if (inter)
        {
            pInter = V0 + alpha*V10 + beta*V20;
            scalar s = (pInter - origin[originI])[i0];

            if ((s >= 0) && (s <= maxLength))
            {
                return true;
            }
        }
    }
    return false;
}



// Intersect triangle with bounding box. Return true if
// any of the faces of bb intersect triangle.
// Note: so returns false if triangle inside bb.
bool Foam::triangleFuncs::intersectBb
(
    const point& p0,
    const point& p1,
    const point& p2,
    const treeBoundBox& cubeBb
)
{
    const vector p10 = p1 - p0;
    const vector p20 = p2 - p0;

    // cubeBb points; counted as if cell with vertex0 at cubeBb.min().
    const point& min = cubeBb.min();
    const point& max = cubeBb.max();

    const point& cube0 = min;
    const point  cube1(min.x(), min.y(), max.z());
    const point  cube2(max.x(), min.y(), max.z());
    const point  cube3(max.x(), min.y(), min.z());

    const point  cube4(min.x(), max.y(), min.z());
    const point  cube5(min.x(), max.y(), max.z());
    const point  cube7(max.x(), max.y(), min.z());

    //
    // Intersect all 12 edges of cube with triangle
    //

    point pInter;
    pointField origin(4);
    // edges in x direction
    origin[0] = cube0;
    origin[1] = cube1;
    origin[2] = cube5;
    origin[3] = cube4;

    scalar maxSx = max.x() - min.x();

    if (intersectAxesBundle(p0, p10, p20, 0, origin, maxSx, pInter))
    {
        return true;
    }

    // edges in y direction
    origin[0] = cube0;
    origin[1] = cube1;
    origin[2] = cube2;
    origin[3] = cube3;

    scalar maxSy = max.y() - min.y();

    if (intersectAxesBundle(p0, p10, p20, 1, origin, maxSy, pInter))
    {
        return true;
    }

    // edges in z direction
    origin[0] = cube0;
    origin[1] = cube3;
    origin[2] = cube7;
    origin[3] = cube4;

    scalar maxSz = max.z() - min.z();

    if (intersectAxesBundle(p0, p10, p20, 2, origin, maxSz, pInter))
    {
        return true;
    }


    // Intersect triangle edges with bounding box
    if (cubeBb.intersects(p0, p1, pInter))
    {
        return true;
    }
    if (cubeBb.intersects(p1, p2, pInter))
    {
        return true;
    }
    if (cubeBb.intersects(p2, p0, pInter))
    {
        return true;
    }

    return false;
}


bool Foam::triangleFuncs::intersect
(
    const point& va0,
    const point& va10,
    const point& va20,

    const point& base,
    const point& normal,

    point& pInter0,
    point& pInter1
)
{
    // Get triangle normal
    vector na = va10 ^ va20;
    scalar magArea = mag(na);
    na/magArea;

    if (mag(na & normal) > (1 - SMALL))
    {
        // Parallel
        return false;
    }

    const point va1 = va0 + va10;
    const point va2 = va0 + va20;

    // Find the triangle point on the other side.
    scalar sign0 = (va0 - base) & normal;
    scalar sign1 = (va1 - base) & normal;
    scalar sign2 = (va2 - base) & normal;

    label oppositeVertex = -1;

    if (sign0 < 0)
    {
        if (sign1 < 0)
        {
            if (sign2 < 0)
            {
                // All on same side of plane
                return false;
            }
            else    // sign2 >= 0
            {
                // 2 on opposite side.
                oppositeVertex = 2;
            }
        }
        else    // sign1 >= 0
        {
            if (sign2 < 0)
            {
                // 1 on opposite side.
                oppositeVertex = 1;
            }
            else
            {
                // 0 on opposite side.
                oppositeVertex = 0;
            }
        }
    }
    else    // sign0 >= 0
    {
        if (sign1 < 0)
        {
            if (sign2 < 0)
            {
                // 0 on opposite side.
                oppositeVertex = 0;
            }
            else    // sign2 >= 0
            {
                // 1 on opposite side.
                oppositeVertex = 1;
            }
        }
        else    // sign1 >= 0
        {
            if (sign2 < 0)
            {
                // 2 on opposite side.
                oppositeVertex = 2;
            }
            else    // sign2 >= 0
            {
                // All on same side of plane
                return false;
            }
        }
    }

    scalar tol = SMALL*Foam::sqrt(magArea);

    if (oppositeVertex == 0)
    {
        // 0 on opposite side. Cut edges 01 and 02
        setIntersection(va0, sign0, va1, sign1, tol, pInter0);
        setIntersection(va0, sign0, va2, sign2, tol, pInter1);
    }
    else if (oppositeVertex == 1)
    {
        // 1 on opposite side. Cut edges 10 and 12
        setIntersection(va1, sign1, va0, sign0, tol, pInter0);
        setIntersection(va1, sign1, va2, sign2, tol, pInter1);
    }
    else // oppositeVertex == 2
    {
        // 2 on opposite side. Cut edges 20 and 21
        setIntersection(va2, sign2, va0, sign0, tol, pInter0);
        setIntersection(va2, sign2, va1, sign1, tol, pInter1);
    }

    return true;
}


bool Foam::triangleFuncs::intersect
(
    const point& va0,
    const point& va10,
    const point& va20,

    const point& vb0,
    const point& vb10,
    const point& vb20,

    point& pInter0,
    point& pInter1
)
{
    // Get triangle normals
    vector na = va10 ^ va20;
    na/mag(na);

    vector nb = vb10 ^ vb20;
    nb/mag(nb);

    // Calculate intersection of triangle a with plane of b
    point planeB0;
    point planeB1;
    if (!intersect(va0, va10, va20, vb0, nb, planeB0, planeB1))
    {
        return false;
    }

    //       ,,  triangle b with plane of a
    point planeA0;
    point planeA1;
    if (!intersect(vb0, vb10, vb20, va0, na, planeA0, planeA1))
    {
        return false;
    }

    // Now check if intersections overlap (w.r.t. intersection of the two
    // planes)

    vector intersection(na ^ nb);

    scalar coordB0 = planeB0 & intersection;
    scalar coordB1 = planeB1 & intersection;

    scalar coordA0 = planeA0 & intersection;
    scalar coordA1 = planeA1 & intersection;

    // Put information in indexable form for sorting.
    List<point*> pts(4);
    boolList isFromB(4);
    SortableList<scalar> sortCoords(4);

    pts[0] = &planeB0;
    isFromB[0] = true;
    sortCoords[0] = coordB0;

    pts[1] = &planeB1;
    isFromB[1] = true;
    sortCoords[1] = coordB1;

    pts[2] = &planeA0;
    isFromB[2] = false;
    sortCoords[2] = coordA0;

    pts[3] = &planeA1;
    isFromB[3] = false;
    sortCoords[3] = coordA1;

    sortCoords.sort();

    const labelList& indices = sortCoords.indices();

    if (isFromB[indices[0]] == isFromB[indices[1]])
    {
        // Entry 0 and 1 are of same region (both a or both b). Hence that
        // region does not overlap.
        return false;
    }
    else
    {
        // Different type. Start of overlap at indices[1], end at indices[2]
        pInter0 = *pts[indices[1]];
        pInter1 = *pts[indices[2]];

        return true;
    }
}


bool Foam::triangleFuncs::classify
(
    const point& baseVertex,
    const vector& E0,
    const vector& E1,
    const vector& n,
    const point& pInter,
    const scalar tol,
    label& nearType,
    label& nearLabel
)
{
    // Get largest component of normal
    scalar magX = Foam::mag(n.x());
    scalar magY = Foam::mag(n.y());
    scalar magZ = Foam::mag(n.z());

    label i0 = -1;
    if ((magX >= magY) && (magX >= magZ))
    {
        i0 = 0;
    }
    else if ((magY >= magX) && (magY >= magZ))
    {
        i0 = 1;
    }
    else
    {
        i0 = 2;
    }

    // Get other components
    label i1 = (i0 + 1) % 3;
    label i2 = (i1 + 1) % 3;

    scalar u1 = E0[i1];
    scalar v1 = E0[i2];

    scalar u2 = E1[i1];
    scalar v2 = E1[i2];

    scalar det = v2*u1 - u2*v1;

    scalar u0 = pInter[i1] - baseVertex[i1];
    scalar v0 = pInter[i2] - baseVertex[i2];

    scalar alpha = 0;
    scalar beta = 0;

    bool hit = false;

    if (Foam::mag(u1) < ROOTVSMALL)
    {
        beta = u0/u2;

        alpha = (v0 - beta*v2)/v1;

        hit = ((alpha >= 0) && ((alpha + beta) <= 1));
    }
    else
    {
        beta = (v0*u1 - u0*v1)/det;

        alpha = (u0 - beta*u2)/u1;

        hit = ((alpha >= 0) && ((alpha + beta) <= 1));
    }

    //
    // Now alpha, beta are the coordinates in the triangle local coordinate
    // system E0, E1
    //

    nearType = NONE;
    nearLabel = -1;

    if (mag(alpha+beta-1) < tol)
    {
        // On edge between vert 1-2 (=edge 1)

        if (mag(alpha) < tol)
        {
            nearType = POINT;
            nearLabel = 2;
        }
        else if (mag(beta) < tol)
        {
            nearType = POINT;
            nearLabel = 1;
        }
        else if ((alpha >= 0) && (alpha <= 1) && (beta >= 0) && (beta <= 1))
        {
            nearType = EDGE;
            nearLabel = 1;
        }
    }
    else if (mag(alpha) < tol)
    {
        // On edge between vert 2-0 (=edge 2)

        if (mag(beta) < tol)
        {
            nearType = POINT;
            nearLabel = 0;
        }
        else if (mag(beta-1) < tol)
        {
            nearType = POINT;
            nearLabel = 2;
        }
        else if ((beta >= 0) && (beta <= 1))
        {
            nearType = EDGE;
            nearLabel = 2;
        }
    }
    else if (mag(beta) < tol)
    {
        // On edge between vert 0-1 (= edge 0)

        if (mag(alpha) < tol)
        {
            nearType = POINT;
            nearLabel = 0;
        }
        else if (mag(alpha-1) < tol)
        {
            nearType = POINT;
            nearLabel = 1;
        }
        else if ((alpha >= 0) && (alpha <= 1))
        {
            nearType = EDGE;
            nearLabel = 0;
        }
    }

    return hit;
}


// ************************************************************************* //
