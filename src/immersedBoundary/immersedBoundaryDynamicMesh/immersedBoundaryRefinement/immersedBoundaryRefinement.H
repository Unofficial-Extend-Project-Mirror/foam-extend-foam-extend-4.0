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
