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

Class
    Foam::triSurfaceDistance

Description
    Distance class to triangular surface

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceDistance_H
#define triSurfaceDistance_H

#include "scalarField.H"
#include "triSurfaceMesh.H"
#include "triSurfaceSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class triSurfaceDistance Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceDistance
{
    // Private data

        //- Reference to triangular surface search algorithm
        const triSurfaceSearch& tss_;

        //- Search span
        const vector span_;

        //- Internal volume flag: distance is positive inside the surface
        const bool internal_;

        //- Iterate distance
        const bool iterate_;


    // Private member functions

        //- Return immersed boundary cell span
        vector cellSpan(const label cellID) const;


public:

    // Constructors

        //- Construct from components
        triSurfaceDistance
        (
            const triSurfaceSearch& tss,
            const vector& span,
            const bool internal,
            const bool iterate
        )
        :
            tss_(tss),
            span_(span),
            internal_(internal),
            iterate_(iterate)
        {}


    //- Destructor - default


    // Member Functions

        // Iterate distance
        inline bool iterateDistance() const
        {
            return iterate_;
        }

        // Return distance to triangular surface
        scalar distance(const point& p) const;

        // Return distance to triangular surface
        tmp<scalarField> distance(const vectorField& p) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
