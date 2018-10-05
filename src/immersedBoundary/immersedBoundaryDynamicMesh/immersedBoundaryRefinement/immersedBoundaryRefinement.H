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
    Foam::immersedBoundaryRefinement

Description
    Selection of refinement cells based on distance from immersed boundary
SourceFiles
    immersedBoundaryRefinement.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryRefinement_H
#define immersedBoundaryRefinement_H

#include "refinementSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class immersedBoundaryRefinement Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryRefinement
:
    public refinementSelection
{
    // Private data

        //- Refinement distance to immersed boundary
        scalar refinementDistance_;

        //- Unrefinement distance to immersed boundary
        scalar unrefinementDistance_;

        //- Refinement distance to immersed boundary inside of the body
        scalar internalRefinementDistance_;

        //- Unefinement distance to immersed boundary inside of the body
        scalar internalUnrefinementDistance_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        immersedBoundaryRefinement(const immersedBoundaryRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const immersedBoundaryRefinement&);


public:

    //- Runtime type information
    TypeName("immersedBoundaryRefinement");


    // Constructors

        //- Construct from components
        immersedBoundaryRefinement
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~immersedBoundaryRefinement();


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
