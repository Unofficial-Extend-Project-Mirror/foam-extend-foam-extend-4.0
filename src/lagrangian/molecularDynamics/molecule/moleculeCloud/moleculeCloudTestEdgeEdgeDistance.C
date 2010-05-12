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

\*----------------------------------------------------------------------------*/

#include "moleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::moleculeCloud::testEdgeEdgeDistance
(
    const edge& eI,
    const edge& eJ
) const
{
    const vector& eJs(mesh_.points()[eJ.start()]);
    const vector& eJe(mesh_.points()[eJ.end()]);

    return testEdgeEdgeDistance(eI, eJs, eJe);
}

bool Foam::moleculeCloud::testEdgeEdgeDistance
(
    const edge& eI,
    const vector& eJs,
    const vector& eJe
) const
{
    vector a(eI.vec(mesh_.points()));
    vector b(eJe - eJs);

    const vector& eIs(mesh_.points()[eI.start()]);

    vector c(eJs - eIs);

    vector crossab = a ^ b;
    scalar magCrossSqr = magSqr(crossab);

    if (magCrossSqr > VSMALL)
    {
        // If the edges are parallel then a point-face
        // search will pick them up

        scalar s = ((c ^ b) & crossab)/magCrossSqr;
        scalar t = ((c ^ a) & crossab)/magCrossSqr;

        // Check for end points outside of range 0..1
        // If the closest point is outside this range
        // a point-face search will have found it.

        return
        (
            s >= 0
         && s <= 1
         && t >= 0
         && t <= 1
         && magSqr(eIs + a*s - eJs - b*t) <= pairPotentials_.rCutMaxSqr()
        );
    }

    return false;
}


// ************************************************************************* //
