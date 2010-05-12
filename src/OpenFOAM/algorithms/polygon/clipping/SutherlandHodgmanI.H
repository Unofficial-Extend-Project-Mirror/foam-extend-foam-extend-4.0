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

Author
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline bool Foam::SutherlandHodgman::isPointVisibleFromCuttingSegment
(
    const point2D& p
) const
{
    // We assume a right-hand rule system
    //
    //                      / P3
    //                     /
    //                    /
    //                   /
    //                  /
    //                 /
    //             P1 ---------------- P2
    //
    // if CrossProduct = (x2 − x1)(y3 − y1) − (y2 − y1)(x3 − x1) >= 0 ,
    // so sin(teta) >= 0, we consider the point p3 visible
    //
    // Ref: http://en.wikipedia.org/wiki/Cross_product, Section Applications,
    // Computational geometry
    //
    //  P1 (x1, y1)
    //  P2 (x2, y2)
    //  P3 (x3, y3)
    //
    //  x1 = clippingPolygon_[currentClipEdgeP1_][0]
    //  y1 = clippingPolygon_[currentClipEdgeP1_][1]
    //  x2 = clippingPolygon_[currentClipEdgeP2_][0]
    //  y2 = clippingPolygon_[currentClipEdgeP2_][1]
    //  x3 = p[0]
    //  y3 = p[1]
    
    // return ((x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1) >= 0);
    return
    (
        (
            clippingPolygon_[currentClipEdgeP2_][0]
          - clippingPolygon_[currentClipEdgeP1_][0]
        )*(p[1] - clippingPolygon_[currentClipEdgeP1_][1])
      - (
          clippingPolygon_[currentClipEdgeP2_][1]
        - clippingPolygon_[currentClipEdgeP1_][1]
        )*(p[0] - clippingPolygon_[currentClipEdgeP1_][0])
     >= 0
    );
}


inline bool Foam::SutherlandHodgman::clipSegment
(
    const point2D& p1,
    const point2D& p2,
    point2D& intersectPt
) const
{
    return lineSegmentIntersection
    (
        p1,
        p2,
        clippingPolygon_[currentClipEdgeP1_],
        clippingPolygon_[currentClipEdgeP2_],
        intersectPt
    );
}


// ************************************************************************* //
