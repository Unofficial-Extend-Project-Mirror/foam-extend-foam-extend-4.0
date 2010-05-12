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

\*---------------------------------------------------------------------------*/

#include "treeBoundBox.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::treeBoundBox Foam::treeBoundBox::greatBox
(
    vector(-GREAT, -GREAT, -GREAT),
    vector(GREAT, GREAT, GREAT)
);


//! @cond - skip documentation : local scope only
const Foam::label facesArray[6][4] =
{
    {0, 4, 6, 2}, // left
    {1, 3, 7, 5}, // right
    {0, 1, 5, 4}, // bottom
    {2, 6, 7, 3}, // top
    {0, 2, 3, 1}, // back
    {4, 5, 7, 6}  // front
};
//! @endcond


const Foam::faceList Foam::treeBoundBox::faces
(
    initListList<face, label, 6, 4>(facesArray)
);


//! @cond - skip documentation : local scope only
const Foam::label edgesArray[12][2] =
{
    {0, 1}, // 0
    {1, 3},
    {2, 3}, // 2
    {0, 2},
    {4, 5}, // 4
    {5, 7},
    {6, 7}, // 6
    {4, 6},
    {0, 4}, // 8
    {1, 5},
    {3, 7}, // 10
    {2, 6}
};
//! @endcond


const Foam::edgeList Foam::treeBoundBox::edges
(
    initListList<edge, label, 12, 2>(edgesArray)
);


const Foam::FixedList<Foam::vector, 6> Foam::treeBoundBox::faceNormals
(
    calcFaceNormals()
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::FixedList<Foam::vector, 6> Foam::treeBoundBox::calcFaceNormals()
{
    FixedList<vector, 6> normals;
    normals[LEFT]   = vector(-1,  0,  0);
    normals[RIGHT]  = vector( 1,  0,  0);
    normals[BOTTOM] = vector( 0, -1,  0);
    normals[TOP]    = vector( 0,  1,  0);
    normals[BACK]   = vector( 0,  0, -1);
    normals[FRONT]  = vector( 0,  0,  1);
    return normals;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct as the bounding box of the given pointField
Foam::treeBoundBox::treeBoundBox(const UList<point>& points)
:
    boundBox()
{
    if (points.size() == 0)
    {
        WarningIn("treeBoundBox::treeBoundBox(const UList<point>&)")
            << "cannot find bounding box for zero sized pointField"
            << "returning zero" << endl;

        return;
    }

    min() = points[0];
    max() = points[0];

    for (label i = 1; i < points.size(); i++)
    {
        min() = ::Foam::min(min(), points[i]);
        max() = ::Foam::max(max(), points[i]);
    }
}


// Construct as the bounding box of the given pointField
Foam::treeBoundBox::treeBoundBox
(
    const UList<point>& points,
    const labelList& meshPoints
)
:
    boundBox()
{
    if (meshPoints.size() == 0)
    {
        WarningIn
        (
            "treeBoundBox::treeBoundBox(const UList<point>&, const labelList)"
        )   << "cannot find bounding box for zero sized pointField"
            << "returning zero" << endl;

        return;
    }

    min() = points[meshPoints[0]];
    max() = points[meshPoints[0]];

    for (label i = 1; i < meshPoints.size(); i++)
    {
        min() = ::Foam::min(min(), points[meshPoints[i]]);
        max() = ::Foam::max(max(), points[meshPoints[i]]);
    }
}


// Construct from Istream
Foam::treeBoundBox::treeBoundBox(Istream& is)
:
    boundBox(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointField Foam::treeBoundBox::points() const
{
    pointField points(8);
    forAll(points, octant)
    {
        points[octant] = corner(octant);

    }
    return points;
}


Foam::treeBoundBox Foam::treeBoundBox::subBbox(const direction octant) const
{
    if (octant > 7)
    {
        FatalErrorIn
        (
            "treeBoundBox::subCube(const direction)"
        )   << "octant should be [0..7]"
            << abort(FatalError);
    }

    scalar leftx, lefty, leftz;
    scalar rightx, righty, rightz;

    scalar midx=0.5*(min().x() + max().x());
    scalar midy=0.5*(min().y() + max().y());
    scalar midz=0.5*(min().z() + max().z());

    // X half
    if (octant & treeBoundBox::RIGHTHALF)
    {
        leftx = midx;
        rightx = max().x();
    }
    else
    {
        leftx = min().x();
        rightx = midx;
    }

    // Y half
    if (octant & treeBoundBox::TOPHALF)
    {
        lefty = midy;
        righty = max().y();
    }
    else
    {
        lefty = min().y();
        righty = midy;
    }

    // Z half
    if (octant & treeBoundBox::FRONTHALF)
    {
        leftz = midz;
        rightz = max().z();
    }
    else
    {
        leftz = min().z();
        rightz = midz;
    }

    return treeBoundBox
    (
        point(leftx, lefty, leftz),
        point(rightx, righty, rightz)
    );
}


// Octant to bounding box using permutation only.
Foam::treeBoundBox Foam::treeBoundBox::subBbox
(
    const point& mid,
    const direction octant
) const
{
    if (octant > 7)
    {
        FatalErrorIn
        (
            "treeBoundBox::subCube(const point&, const direction)"
        )   << "octant should be [0..7]"
            << abort(FatalError);
    }

    treeBoundBox subBb;
    point& subMin = subBb.min();
    point& subMax = subBb.max();

    if (octant & treeBoundBox::RIGHTHALF)
    {
        subMin.x() = mid.x();
        subMax.x() = max().x();
    }
    else
    {
        subMin.x() = min().x();
        subMax.x() = mid.x();
    }
    if (octant & treeBoundBox::TOPHALF)
    {
        subMin.y() = mid.y();
        subMax.y() = max().y();
    }
    else
    {
        subMin.y() = min().y();
        subMax.y() = mid.y();
    }
    if (octant & treeBoundBox::FRONTHALF)
    {
        subMin.z() = mid.z();
        subMax.z() = max().z();
    }
    else
    {
        subMin.z() = min().z();
        subMax.z() = mid.z();
    }

    return subBb;
}

bool Foam::treeBoundBox::overlaps
(
    const point& centre,
    const scalar radiusSqr
) const
{
    // Find out where centre is in relation to bb.
    // Find nearest point on bb.
    scalar distSqr = 0;

    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        scalar d0 = min()[dir] - centre[dir];
        scalar d1 = max()[dir] - centre[dir];

        if ((d0 > 0) != (d1 > 0))
        {
            // centre inside both extrema. This component does not add any
            // distance.
        }
        else if (Foam::mag(d0) < Foam::mag(d1))
        {
            distSqr += d0*d0;
        }
        else
        {
            distSqr += d1*d1;
        }

        if (distSqr > radiusSqr)
        {
            return false;
        }
    }

    return true;
}
// line intersection. Returns true if line (start to end) inside
// bb or intersects bb. Sets pt to intersection.
//
// Sutherlands algorithm:
//   loop
//     - start = intersection of line with one of the planes bounding
//       the bounding box
//     - stop if start inside bb (return true)
//     - stop if start and end in same 'half' (e.g. both above bb)
//       (return false)
//
// Uses posBits to efficiently determine 'half' in which start and end
// point are.
//
// Note:
//   - sets coordinate to exact position: e.g. pt.x() = min().x()
//     since plane intersect routine might have truncation error.
//     This makes sure that posBits tests 'inside'
bool Foam::treeBoundBox::intersects
(
    const point& start,
    const point& end,
    point& pt
) const
{
    vector vec(end - start);

    pt = start;

    const direction endBits = posBits(end);

    while(true)
    {
        direction ptBits = posBits(pt);

        if (ptBits == 0)
        {
            // pt inside bb
            return true;
        }

        if ((ptBits & endBits) != 0)
        {
            // pt and end in same block outside of bb
            return false;
        }

        if (ptBits & LEFTBIT)
        {
            // Intersect with plane V=min, n=-1,0,0
            if (Foam::mag(vec.x()) > VSMALL)
            {
                scalar s = (min().x() - pt.x())/vec.x();
                pt.x() = min().x();
                pt.y() = pt.y() + vec.y()*s;
                pt.z() = pt.z() + vec.z()*s;
            }
        }
        if (ptBits & RIGHTBIT)
        {
            // Intersect with plane V=max, n=1,0,0
            if (Foam::mag(vec.x()) > VSMALL)
            {
                scalar s = (max().x() - pt.x())/vec.x();
                pt.x() = max().x();
                pt.y() = pt.y() + vec.y()*s;
                pt.z() = pt.z() + vec.z()*s;
            }
        }

        if (ptBits & BOTTOMBIT)
        {
            // Intersect with plane V=min, n=0,-1,0
            if (Foam::mag(vec.y()) > VSMALL)
            {
                scalar s = (min().y() - pt.y())/vec.y();
                pt.x() = pt.x() + vec.x()*s;
                pt.y() = min().y();
                pt.z() = pt.z() + vec.z()*s;
            }
        }
        if (ptBits & TOPBIT)
        {
            // Intersect with plane V=max, n=0,1,0
            if (Foam::mag(vec.y()) > VSMALL)
            {
                scalar s = (max().y() - pt.y())/vec.y();
                pt.x() = pt.x() + vec.x()*s;
                pt.y() = max().y();
                pt.z() = pt.z() + vec.z()*s;
            }
        }

        if (ptBits & BACKBIT)
        {
            // Intersect with plane V=min, n=0,0,-1
            if (Foam::mag(vec.z()) > VSMALL)
            {
                scalar s = (min().z() - pt.z())/vec.z();
                pt.x() = pt.x() + vec.x()*s;
                pt.y() = pt.y() + vec.y()*s;
                pt.z() = min().z();
            }
        }
        if (ptBits & FRONTBIT)
        {
            // Intersect with plane V=max, n=0,0,1
            if (Foam::mag(vec.z()) > VSMALL)
            {
                scalar s = (max().z() - pt.z())/vec.z();
                pt.x() = pt.x() + vec.x()*s;
                pt.y() = pt.y() + vec.y()*s;
                pt.z() = max().z();
            }
        }
    }
}


// this.bb fully contains bb
bool Foam::treeBoundBox::contains(const treeBoundBox& bb) const
{
    return contains(bb.min()) && contains(bb.max());
}


bool Foam::treeBoundBox::containsNarrow(const point& sample) const
{
    return
    (
        (sample.x() > min().x()) &&
        (sample.y() > min().y()) &&
        (sample.z() > min().z()) &&
        (sample.x() < max().x()) &&
        (sample.y() < max().y()) &&
        (sample.z() < max().z())
    );
}

bool Foam::treeBoundBox::contains(const vector& dir, const point& sample) const
{
    //
    // Compare all components against min and max of bb
    //

    for (direction cmpt=0; cmpt<3; cmpt++)
    {
        if (sample[cmpt] < min()[cmpt])
        {
            return false;
        }
        else if (sample[cmpt] == min()[cmpt])
        {
            // On edge. Outside if direction points outwards.
            if (dir[cmpt] < 0)
            {
                return false;
            }
        }

        if (sample[cmpt] > max()[cmpt])
        {
            return false;
        }
        else if (sample[cmpt] == max()[cmpt])
        {
            // On edge. Outside if direction points outwards.
            if (dir[cmpt] > 0)
            {
                return false;
            }
        }
    }

    // All components inside bb
    return true;
}


// Code position of point relative to box
Foam::direction Foam::treeBoundBox::posBits(const point& pt) const
{
    direction posBits = 0;

    if (pt.x() < min().x())
    {
        posBits |= LEFTBIT;
    }
    if (pt.x() > max().x())
    {
        posBits |= RIGHTBIT;
    }

    if (pt.y() < min().y())
    {
        posBits |= BOTTOMBIT;
    }
    if (pt.y() > max().y())
    {
        posBits |= TOPBIT;
    }

    if (pt.z() < min().z())
    {
        posBits |= BACKBIT;
    }
    if (pt.z() > max().z())
    {
        posBits |= FRONTBIT;
    }
    return posBits;
}


// nearest and furthest corner coordinate.
// !names of treeBoundBox::min() and treeBoundBox::max() are confusing!
void Foam::treeBoundBox::calcExtremities
(
    const point& sample,
    point& nearest,
    point& furthest
) const
{
    scalar nearX, nearY, nearZ;
    scalar farX, farY, farZ;

    if (Foam::mag(min().x() - sample.x()) < Foam::mag(max().x() - sample.x()))
    {
        nearX = min().x();
        farX = max().x();
    }
    else
    {
        nearX = max().x();
        farX = min().x();
    }

    if (Foam::mag(min().y() - sample.y()) < Foam::mag(max().y() - sample.y()))
    {
        nearY = min().y();
        farY = max().y();
    }
    else
    {
        nearY = max().y();
        farY = min().y();
    }

    if (Foam::mag(min().z() - sample.z()) < Foam::mag(max().z() - sample.z()))
    {
        nearZ = min().z();
        farZ = max().z();
    }
    else
    {
        nearZ = max().z();
        farZ = min().z();
    }

    nearest = point(nearX, nearY, nearZ);
    furthest = point(farX, farY, farZ);
}


Foam::scalar Foam::treeBoundBox::maxDist(const point& sample) const
{
    point near, far;
    calcExtremities(sample, near, far);

    return Foam::mag(far - sample);
}


// Distance comparator
// Compare all vertices of bounding box against all of other bounding
// box to see if all vertices of one are nearer
Foam::label Foam::treeBoundBox::distanceCmp
(
    const point& sample,
    const treeBoundBox& other
) const
{
    //
    // Distance sample <-> nearest and furthest away vertex of this
    //

    point nearThis, farThis;

    // get nearest and furthest away vertex
    calcExtremities(sample, nearThis, farThis);

    const scalar minDistThis =
        sqr(nearThis.x() - sample.x())
     +  sqr(nearThis.y() - sample.y())
     +  sqr(nearThis.z() - sample.z());
    const scalar maxDistThis =
        sqr(farThis.x() - sample.x())
     +  sqr(farThis.y() - sample.y())
     +  sqr(farThis.z() - sample.z());

    //
    // Distance sample <-> other
    //

    point nearOther, farOther;

    // get nearest and furthest away vertex
    other.calcExtremities(sample, nearOther, farOther);

    const scalar minDistOther =
        sqr(nearOther.x() - sample.x())
     +  sqr(nearOther.y() - sample.y())
     +  sqr(nearOther.z() - sample.z());
    const scalar maxDistOther =
        sqr(farOther.x() - sample.x())
     +  sqr(farOther.y() - sample.y())
     +  sqr(farOther.z() - sample.z());

    //
    // Categorize
    //
    if (maxDistThis < minDistOther)
    {
        // All vertices of this are nearer to sample than any vertex of other
        return -1;
    }
    else if (minDistThis > maxDistOther)
    {
        // All vertices of this are further from sample than any vertex of other
        return 1;
    }
    else
    {
        // Mixed bag
        return 0;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool Foam::operator==(const treeBoundBox& a, const treeBoundBox& b)
{
    return (a.min() == b.min()) && (a.max() == b.max());
}


bool Foam::operator!=(const treeBoundBox& a, const treeBoundBox& b)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * IOstream Operator  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, treeBoundBox& bb)
{
    is >> bb.min() >> bb.max();
    return is;
}


// ************************************************************************* //
