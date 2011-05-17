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

Implemented by
    Sandeep Menon
    University of Massachusetts Amherst

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Compute clip-planes
inline void tetIntersection::computeClipPlanes()
{
    // Define edge vectors
    vector edge10 = clipTet_[1] - clipTet_[0];
    vector edge20 = clipTet_[2] - clipTet_[0];
    vector edge30 = clipTet_[3] - clipTet_[0];
    vector edge21 = clipTet_[2] - clipTet_[1];
    vector edge31 = clipTet_[3] - clipTet_[1];

    // Cross-products
    clipPlanes_[0].first() = (edge20 ^ edge10);
    clipPlanes_[1].first() = (edge10 ^ edge30);
    clipPlanes_[2].first() = (edge30 ^ edge20);
    clipPlanes_[3].first() = (edge21 ^ edge31);

    // Normalize
    clipPlanes_[0].first() /= mag(clipPlanes_[0].first()) + VSMALL;
    clipPlanes_[1].first() /= mag(clipPlanes_[1].first()) + VSMALL;
    clipPlanes_[2].first() /= mag(clipPlanes_[2].first()) + VSMALL;
    clipPlanes_[3].first() /= mag(clipPlanes_[3].first()) + VSMALL;

    // Compute magnitude of clipping tetrahedron
    clipTetMag_ = (1.0 / 6.0) * (edge10 & clipPlanes_[3].first());

    if (clipTetMag_ < 0.0)
    {
        // Reverse normal directions
        clipPlanes_[0].first() = -clipPlanes_[0].first();
        clipPlanes_[1].first() = -clipPlanes_[1].first();
        clipPlanes_[2].first() = -clipPlanes_[2].first();
        clipPlanes_[3].first() = -clipPlanes_[3].first();

        // Reverse sign
        clipTetMag_ = mag(clipTetMag_);
    }

    // Determine plane constants
    clipPlanes_[0].second() = (clipTet_[0] & clipPlanes_[0].first());
    clipPlanes_[1].second() = (clipTet_[1] & clipPlanes_[1].first());
    clipPlanes_[2].second() = (clipTet_[2] & clipPlanes_[2].first());
    clipPlanes_[3].second() = (clipTet_[3] & clipPlanes_[3].first());
}


// Split and decompose
inline void tetIntersection::splitAndDecompose
(
    const label tetPlaneIndex,
    FixedList<point, 4>& tetra,
    DynamicList<FixedList<point, 4> >& decompTets
) const
{
    FixedList<scalar, 4> C;
    FixedList<vector, 4> tmpTetra;
    FixedList<label, 4> pos, neg, zero;
    label i = 0, nPos = 0, nNeg = 0, nZero = 0;

    // Fetch reference to plane
    const hPlane& tetPlane = clipPlanes_[tetPlaneIndex];

    for (i = 0; i < 4; ++i)
    {
        // Compute distance to plane
        C[i] = (tetra[i] & tetPlane.first()) - tetPlane.second();

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
        decompTets.append(tetra);
        return;
    }

    // Tetrahedron is split by plane.  Determine how it is split and how to
    // decompose the negative-side portion into tetrahedra (6 cases).
    scalar w0, w1, invCDiff;
    vector intp[4];

    if (nPos == 3)
    {
        // +++-
        for (i = 0; i < nPos; ++i)
        {
            invCDiff = (1.0 / (C[pos[i]] - C[neg[0]]));

            w0 = -C[neg[0]] * invCDiff;
            w1 = +C[pos[i]] * invCDiff;

            tetra[pos[i]] = (w0 * tetra[pos[i]]) + (w1 * tetra[neg[0]]);
        }

        decompTets.append(tetra);
    }
    else
    if (nPos == 2)
    {
        if (nNeg == 2)
        {
            // ++--
            for (i = 0; i < nPos; ++i)
            {
                invCDiff = (1.0 / (C[pos[i]] - C[neg[0]]));

                w0 = -C[neg[0]] * invCDiff;
                w1 = +C[pos[i]] * invCDiff;

                intp[i] = (w0 * tetra[pos[i]]) + (w1 * tetra[neg[0]]);
            }

            for (i = 0; i < nNeg; ++i)
            {
                invCDiff = (1.0 / (C[pos[i]] - C[neg[1]]));

                w0 = -C[neg[1]] * invCDiff;
                w1 = +C[pos[i]] * invCDiff;

                intp[i+2] = (w0 * tetra[pos[i]]) + (w1 * tetra[neg[1]]);
            }

            tetra[pos[0]] = intp[2];
            tetra[pos[1]] = intp[1];

            decompTets.append(tetra);

            tmpTetra[0] = tetra[neg[1]];
            tmpTetra[1] = intp[3];
            tmpTetra[2] = intp[2];
            tmpTetra[3] = intp[1];

            decompTets.append(tmpTetra);

            tmpTetra[0] = tetra[neg[0]];
            tmpTetra[1] = intp[0];
            tmpTetra[2] = intp[1];
            tmpTetra[3] = intp[2];

            decompTets.append(tmpTetra);
        }
        else
        {
            // ++-0
            for (i = 0; i < nPos; ++i)
            {
                invCDiff = (1.0 / (C[pos[i]] - C[neg[0]]));

                w0 = -C[neg[0]] * invCDiff;
                w1 = +C[pos[i]] * invCDiff;

                tetra[pos[i]] = (w0 * tetra[pos[i]]) + (w1 * tetra[neg[0]]);
            }

            decompTets.append(tetra);
        }
    }
    else
    if (nPos == 1)
    {
        if (nNeg == 3)
        {
            // +---
            for (i = 0; i < nNeg; ++i)
            {
                invCDiff = (1.0 / (C[pos[0]] - C[neg[i]]));

                w0 = -C[neg[i]] * invCDiff;
                w1 = +C[pos[0]] * invCDiff;

                intp[i] = (w0 * tetra[pos[0]]) + (w1 * tetra[neg[i]]);
            }

            tetra[pos[0]] = intp[0];

            decompTets.append(tetra);

            tmpTetra[0] = intp[0];
            tmpTetra[1] = tetra[neg[1]];
            tmpTetra[2] = tetra[neg[2]];
            tmpTetra[3] = intp[1];

            decompTets.append(tmpTetra);

            tmpTetra[0] = tetra[neg[2]];
            tmpTetra[1] = intp[1];
            tmpTetra[2] = intp[2];
            tmpTetra[3] = intp[0];

            decompTets.append(tmpTetra);
        }
        else
        if (nNeg == 2)
        {
            // +--0
            for (i = 0; i < nNeg; ++i)
            {
                invCDiff = (1.0 / (C[pos[0]] - C[neg[i]]));

                w0 = -C[neg[i]] * invCDiff;
                w1 = +C[pos[0]] * invCDiff;

                intp[i] = (w0 * tetra[pos[0]]) + (w1 * tetra[neg[i]]);
            }

            tetra[pos[0]] = intp[0];

            decompTets.append(tetra);

            tmpTetra[0] = intp[1];
            tmpTetra[1] = tetra[zero[0]];
            tmpTetra[2] = tetra[neg[1]];
            tmpTetra[3] = intp[0];

            decompTets.append(tmpTetra);
        }
        else
        {
            // +-00
            invCDiff = (1.0 / (C[pos[0]] - C[neg[0]]));

            w0 = -C[neg[0]] * invCDiff;
            w1 = +C[pos[0]] * invCDiff;

            tetra[pos[0]] = (w0 * tetra[pos[0]]) + (w1 * tetra[neg[0]]);

            decompTets.append(tetra);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline tetIntersection::tetIntersection(const FixedList<point, 4>& clipTet)
:
    clipTet_(clipTet),
    clipTetMag_(0.0),
    inside_(10),
    allTets_(10)
{
    // Pre-compute clipping planes
    computeClipPlanes();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline tetIntersection::~tetIntersection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return magnitude of clipping tetrahedron
inline scalar tetIntersection::clipTetMag() const
{
    return clipTetMag_;
}


// Evaluate for intersections
inline bool tetIntersection::evaluate(const FixedList<point, 4>& subjectTet)
{
    // Clear lists
    inside_.clear();
    allTets_.clear();

    // Add initial tetrahedron to list
    allTets_.append(subjectTet);

    // Clip second tetrahedron against planes of clipping tetrahedron
    for (label i = 0; i < 4; i++)
    {
        forAll(allTets_, tetI)
        {
            splitAndDecompose(i, allTets_[tetI], inside_);
        }

        // Prep for next clipping plane
        allTets_ = inside_;
        inside_.clear();
    }

    return (allTets_.size() > 0);
}


// Return intersections
inline const DynamicList<FixedList<point, 4> >&
tetIntersection::getIntersection() const
{
    return allTets_;
}


//- Evaluate and return volume / centroid
inline void tetIntersection::getVolumeAndCentre
(
    scalar& volume,
    vector& centre
) const
{
    volume = 0.0;
    centre = vector::zero;

    forAll(allTets_, tetI)
    {
        const FixedList<point, 4>& t = allTets_[tetI];

        // Calculate volume (no check for orientation)
        scalar tV =
        (
            Foam::mag
            (
                (1.0/6.0) *
                (
                    ((t[1] - t[0]) ^ (t[2] - t[0])) & (t[3] - t[0])
                )
            )
        );

        // Calculate centroid
        vector tC = (0.25 * (t[0] + t[1] + t[2] + t[3]));

        volume += tV;
        centre += (tV * tC);
    }

    centre /= volume + VSMALL;
}


}

// ************************************************************************* //
