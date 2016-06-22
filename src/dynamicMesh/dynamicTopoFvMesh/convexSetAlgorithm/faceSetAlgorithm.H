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
    faceSetAlgorithm

Description
    Class for convexSetAlgorithms pertaining to faces

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    faceSetAlgorithm.C

\*---------------------------------------------------------------------------*/

#ifndef faceSetAlgorithm_H
#define faceSetAlgorithm_H

#include "convexSetAlgorithm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class faceSetAlgorithm Declaration
\*---------------------------------------------------------------------------*/

class faceSetAlgorithm
:
    public convexSetAlgorithm
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        faceSetAlgorithm(const faceSetAlgorithm&);

        //- Disallow default bitwise assignment
        void operator=(const faceSetAlgorithm&);

        //- Write out tris as a VTK
        void writeVTK
        (
            const word&,
            const List<FixedList<point, 3> >&
        ) const;

public:

    //- Constructor

        faceSetAlgorithm
        (
            const polyMesh& mesh,
            const pointField& newPoints,
            const UList<edge>& newEdges,
            const UList<face>& newFaces,
            const UList<cell>& newCells,
            const UList<label>& newOwner,
            const UList<label>& newNeighbour
        );

    //- Destructor

        virtual ~faceSetAlgorithm()
        {}

    //- Virtual member functions

        // Dimensions of the algorithm
        virtual label dimension() const
        {
            return 2;
        }

        // Compute normFactor
        virtual void computeNormFactor(const label index) const;

        // Check whether the bounding box contains the entity
        virtual bool contains(const label index) const;

        // Compute intersections
        virtual bool computeIntersection
        (
            const label newIndex,
            const label oldIndex,
            const label offset,
            bool output
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
