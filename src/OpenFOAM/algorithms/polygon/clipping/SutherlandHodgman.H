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

Class
    SutherlandHodgman

Description
    This is an implementation of the Sutherland Hodgman algorithm:
    Reentrant Polygon Clipping, Sutherland, Hodgman,
    Communications of the ACM, 1974
  
    Wikipedia has a very simple Pseudo-Code example of this algorithm
    http://en.wikipedia.org/wiki/Sutherland-Hodgeman.
  
    The subject polygon will be clipped by the clipping polygon. 

Author
    Martin Beaudoin, Hydro-Quebec, (2008)

SourceFiles
    SutherlandHodgmanI.H
    SutherlandHodgman.C

\*---------------------------------------------------------------------------*/

#ifndef SutherlandHodgman_H
#define SutherlandHodgman_H

#include "DynamicList.H"
#include "point2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class SutherlandHodgman Declaration
\*---------------------------------------------------------------------------*/

class SutherlandHodgman
{
    // Private data

        //- The current subjectPolygon
        List<point2D> subjectPolygon_;

        //- The clipping polygon
        List<point2D> clippingPolygon_;

        //- Clipped edge, polygon 1
        label currentClipEdgeP1_;

        //- Clipped edge, polygon 2
        label currentClipEdgeP2_;

        //- Tolerance for detecting intersecting segment
        scalar intersectSegDistTol_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        SutherlandHodgman(const SutherlandHodgman&);

        //- Disallow default bitwise assignment
        void operator=(const SutherlandHodgman&);


        //- Compute intersection between 2 line segments
        bool lineSegmentIntersection
        (
            const point2D& p1,
            const point2D& p2,
            const point2D& p3,
            const point2D& p4,
            point2D& intersectionPoint
        ) const;

        //- Evaluate if point is visible from current clipping segment
        inline bool isPointVisibleFromCuttingSegment
        (
            const point2D& p
        ) const;

        //- Compute intersection between clipping polygon current segment
        //  and the line formed by the point p1, p2
        inline bool clipSegment
        (
            const point2D& p1,
            const point2D& p2,
            point2D& intersectionPoint
        ) const;


public:

    // Constructors

        //- Construct from component 
        SutherlandHodgman
        (
            const List<point2D>& clippingPolygon,
            const List<point2D>& subjectPolygon,
            const scalar& intersectSegDistTol
        );


    // Destructor - default


    // Member Functions

        // Implements the actual algorithm.
        List<point2D> evaluate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SutherlandHodgmanI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
