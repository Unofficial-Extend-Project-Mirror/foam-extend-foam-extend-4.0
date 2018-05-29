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
    Foam::minCellSizeRefinement

Description
    Selection of refinement cells based on a minimum cell size. Assumes that
    cells are mostly cubic and compares the minimum cell size with V^1/3, where
    V is the cell volume. If the cell size is larger than the specified one,
    cell gets selected for refinement.

SourceFiles
    minCellSizeRefinement.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef minCellSizeRefinement_H
#define minCellSizeRefinement_H

#include "refinementSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class minCellSizeRefinement Declaration
\*---------------------------------------------------------------------------*/

class minCellSizeRefinement
:
    public refinementSelection
{
    // Private data

        //- Minimum cell size
        scalar minDelta_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        minCellSizeRefinement(const minCellSizeRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const minCellSizeRefinement&);


public:

    //- Runtime type information
    TypeName("minCellSizeRefinement");


    // Constructors

        //- Construct from components
        minCellSizeRefinement(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~minCellSizeRefinement();


    // Member Functions

        // Selection of refinement/unrefinement candidates

            //- Return transferable list of cells to refine
            virtual Xfer<labelList> refinementCellCandidates() const;

            //- Return transferable list of split points to unrefine
            virtual Xfer<labelList> unrefinementPointCandidates() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
