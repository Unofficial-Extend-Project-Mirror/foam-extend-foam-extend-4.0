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
    cellSetAlgorithm

Description
    Class for convexSetAlgorithms pertaining to cells

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    cellSetAlgorithm.C

\*---------------------------------------------------------------------------*/

#ifndef cellSetAlgorithm_H
#define cellSetAlgorithm_H

#include "convexSetAlgorithm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class cellSetAlgorithm Declaration
\*---------------------------------------------------------------------------*/

class cellSetAlgorithm
:
    public convexSetAlgorithm
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellSetAlgorithm(const cellSetAlgorithm&);

        //- Disallow default bitwise assignment
        void operator=(const cellSetAlgorithm&);

        //- Write out tets as a VTK
        void writeVTK
        (
            const word&,
            const List<FixedList<point, 4> >&
        ) const;

public:

    //- Constructor

        cellSetAlgorithm
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

        virtual ~cellSetAlgorithm()
        {}

    //- Virtual member functions

        // Dimensions of the algorithm
        virtual label dimension() const
        {
            return 3;
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
