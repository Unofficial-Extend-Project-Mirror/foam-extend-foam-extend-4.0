/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

\*---------------------------------------------------------------------------*/

#include "triSurfaceDistance.H"
#include "triSurfaceTools.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::triSurfaceDistance::distance(const point& p) const
{
    // Find the nearest point on the surface within the span
    const pointIndexHit pih = tss_.nearest(p, span_);

    if (pih.hit())
    {
        // Calculate vector to intersection
        const vector ps = pih.hitPoint() - p;

        // Calculate surface norm for the intersection points
        const vector n =
            triSurfaceTools::surfaceNormal
            (
                tss_.surface(),
                pih.index(),
                pih.hitPoint()
            );

        // Calculate distance to intesection from the inside of the surface
        // Note: "inside" implies negative distance, in order to calculate
        // the "wet" face or cell
        // HJ, 24/Nov/2017
        scalar dist = -sign(ps & n)*mag(ps);

        // Adjust sign if solving for external flows
        if (!internal_)
        {
            dist *= -1;
        }

        return dist;
    }
    else
    {
        FatalErrorIn
        (
            "scalar triSurfaceDistance::distance(const point& p) const"
        )   << "Cannot find nearest distance for point " << p
            << abort(FatalError);

        // Return some distance.  Experimental
        return -mag(pih.rawPoint() - p)*tss_.calcInside(p);
    }
}


Foam::tmp<Foam::scalarField>
Foam::triSurfaceDistance::distance(const vectorField& p) const
{
    tmp<scalarField> tdist(new scalarField(p.size()));
    scalarField& dist = tdist();

    forAll (p, i)
    {
        dist[i] = distance(p[i]);
    }

    return tdist;
}


// ************************************************************************* //
