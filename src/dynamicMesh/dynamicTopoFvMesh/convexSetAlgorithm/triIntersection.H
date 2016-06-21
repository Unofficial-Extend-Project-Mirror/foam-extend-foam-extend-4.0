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

Class
    triIntersection

Description
    Implementation of the triangle intersection algorithm

Implemented by
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    triIntersectionI.H

\*---------------------------------------------------------------------------*/

#ifndef triIntersection_H
#define triIntersection_H

#include "point.H"
#include "label.H"
#include "Tuple2.H"
#include "FixedList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class triIntersection Declaration
\*---------------------------------------------------------------------------*/

class triIntersection
{
    // Private data

        //- Const reference to clipping triangle
        const FixedList<point, 3>& clipTri_;

        //- Clip triangle normal
        vector tNorm_;

        //- Hessian-normal plane definition
        typedef Tuple2<vector, scalar> hPlane;

        FixedList<hPlane, 3> clipPlanes_;

        //- Triangles used as temporaries
        DynamicList<FixedList<point, 3> > inside_;

        //- All intersection tris
        DynamicList<FixedList<point, 3> > allTris_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        triIntersection(const triIntersection&);

        //- Disallow default bitwise assignment
        void operator=(const triIntersection&);

        //- Compute clip-planes
        inline void computeClipPlanes();

        //- Split and decompose
        inline void splitAndDecompose
        (
            const label triPlaneIndex,
            FixedList<point, 3>& tri,
            DynamicList<FixedList<point, 3> >& decompTris
        ) const;

public:

    // Constructors

        //- Construct from components
        inline triIntersection(const FixedList<point, 3>& clipTri);

    // Destructor

        inline ~triIntersection();

    // Member Functions

        //- Evaluate for intersections against input triangle
        inline bool evaluate(const FixedList<point, 3>& subjectTri);

        //- Return intersections
        inline const DynamicList<FixedList<point, 3> >& getIntersection() const;

        //- Evaluate and return area / centroid
        inline void getAreaAndCentre(scalar& area, vector& centre) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triIntersectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
