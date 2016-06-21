/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Implemented by
    Sandeep Menon
    University of Massachusetts Amherst

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Compute clip-planes
inline void triIntersection::computeClipPlanes()
{
    // Face normal
    tNorm_ = (clipTri_[1] - clipTri_[0]) ^ (clipTri_[2] - clipTri_[0]);

    // Normalize
    tNorm_ /= mag(tNorm_) + VSMALL;

    // Edge normals
    clipPlanes_[0].first() = ((clipTri_[1] - clipTri_[0]) ^ tNorm_);
    clipPlanes_[1].first() = ((clipTri_[2] - clipTri_[1]) ^ tNorm_);
    clipPlanes_[2].first() = ((clipTri_[0] - clipTri_[2]) ^ tNorm_);

    // Normalize
    clipPlanes_[0].first() /= mag(clipPlanes_[0].first()) + VSMALL;
    clipPlanes_[1].first() /= mag(clipPlanes_[1].first()) + VSMALL;
    clipPlanes_[2].first() /= mag(clipPlanes_[2].first()) + VSMALL;

    // Determine plane constants
    clipPlanes_[0].second() = (clipTri_[0] & clipPlanes_[0].first());
    clipPlanes_[1].second() = (clipTri_[1] & clipPlanes_[1].first());
    clipPlanes_[2].second() = (clipTri_[2] & clipPlanes_[2].first());
}


// Split and decompose
inline void triIntersection::splitAndDecompose
(
    const label triPlaneIndex,
    FixedList<point, 3>& tri,
    DynamicList<FixedList<point, 3> >& decompTris
) const
{
    FixedList<scalar, 3> C;
    FixedList<vector, 3> tmpTri;
    FixedList<label, 3> pos, neg, zero;
    label i = 0, nPos = 0, nNeg = 0, nZero = 0;

    // Fetch reference to plane
    const hPlane& triPlane = clipPlanes_[triPlaneIndex];

    for (i = 0; i < 3; ++i)
    {
        // Compute distance to plane
        C[i] = (tri[i] & triPlane.first()) - triPlane.second();

        if (C[i] > 0.0)
        {
            pos[nPos++] = i;
        }
        else
        if (C[i] < 0.0)
        {
            neg[nNeg++] = i;
        }
        else
        {
            zero[nZero++] = i;
        }
    }

    if (nNeg == 0)
    {
        return;
    }

    if (nPos == 0)
    {
        decompTris.append(tri);
        return;
    }

    // Triangle is split by plane.  Determine how it is split and how to
    // decompose the negative-side portion into triangles (3 cases).
    scalar w0, w1, invCDiff;
    vector intp[3];

    if (nPos == 2)
    {
        // ++-
        for (i = 0; i < nPos; ++i)
        {
            invCDiff = (1.0 / (C[pos[i]] - C[neg[0]]));

            w0 = -C[neg[0]] * invCDiff;
            w1 = +C[pos[i]] * invCDiff;

            tri[pos[i]] = (w0 * tri[pos[i]]) + (w1 * tri[neg[0]]);
        }

        decompTris.append(tri);
    }
    else
    if (nPos == 1)
    {
        if (nNeg == 2)
        {
            // +--
            for (i = 0; i < nNeg; ++i)
            {
                invCDiff = (1.0 / (C[pos[0]] - C[neg[i]]));

                w0 = -C[neg[i]] * invCDiff;
                w1 = +C[pos[0]] * invCDiff;

                intp[i] = (w0 * tri[pos[0]]) + (w1 * tri[neg[i]]);
            }

            tri[pos[0]] = intp[0];

            decompTris.append(tri);

            tmpTri[0] = intp[0];
            tmpTri[1] = tri[neg[1]];
            tmpTri[2] = intp[1];

            decompTris.append(tmpTri);
        }
        else
        {
            // +-0
            invCDiff = (1.0 / (C[pos[0]] - C[neg[0]]));

            w0 = -C[neg[0]] * invCDiff;
            w1 = +C[pos[0]] * invCDiff;

            tri[pos[0]] = (w0 * tri[pos[0]]) + (w1 * tri[neg[0]]);

            decompTris.append(tri);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline triIntersection::triIntersection(const FixedList<point, 3>& clipTri)
:
    clipTri_(clipTri),
    inside_(10),
    allTris_(10)
{
    // Pre-compute clipping planes
    computeClipPlanes();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline triIntersection::~triIntersection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Evaluate for intersections
inline bool triIntersection::evaluate(const FixedList<point, 3>& subjectTri)
{
    // Clear lists
    inside_.clear();
    allTris_.clear();

    // Project subject triangle to clipping triangle plane
    FixedList<point, 3> sTP(subjectTri);

    for (label i = 0; i < 3; i++)
    {
        // Relative to zeroth face-point
        sTP[i] -= clipTri_[0];

        // Project to face
        sTP[i] = clipTri_[0] + (sTP[i] - (sTP[i] & tNorm_)*tNorm_);
    }

    // Add initial triangle to list
    allTris_.append(sTP);

    // Clip second triangle against planes of clipping triangle
    for (label i = 0; i < 3; i++)
    {
        forAll(allTris_, triI)
        {
            splitAndDecompose(i, allTris_[triI], inside_);
        }

        // Prep for next clipping plane
        allTris_ = inside_;
        inside_.clear();
    }

    return (allTris_.size() > 0);
}


// Return intersections
inline const DynamicList<FixedList<point, 3> >&
triIntersection::getIntersection() const
{
    return allTris_;
}


//- Evaluate and return area / centroid
inline void triIntersection::getAreaAndCentre
(
    scalar& area,
    vector& centre
) const
{
    area = 0.0;
    centre = vector::zero;

    forAll(allTris_, triI)
    {
        const FixedList<point, 3>& t = allTris_[triI];

        // Calculate area (no check for orientation)
        scalar tA = Foam::mag(0.5 * ((t[1] - t[0]) ^ (t[2] - t[0])));

        // Calculate centroid
        vector tC = (1.0 / 3.0) * (t[0] + t[1] + t[2]);

        area += tA;
        centre += (tA * tC);
    }

    centre /= area + VSMALL;
}


}

// ************************************************************************* //
