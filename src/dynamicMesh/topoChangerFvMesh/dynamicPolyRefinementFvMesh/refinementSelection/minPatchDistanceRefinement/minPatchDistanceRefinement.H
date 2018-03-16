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
    Foam::minPatchDistanceRefinement

Description
    Selection of refinement cells based on a minimum distance from a set of
    patches. All cells farther away from the minimum distance get selected.

SourceFiles
    minPatchDistanceRefinement.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef minPatchDistanceRefinement_H
#define minPatchDistanceRefinement_H

#include "refinementSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class minPatchDistanceRefinement Declaration
\*---------------------------------------------------------------------------*/

class minPatchDistanceRefinement
:
    public refinementSelection
{
    // Private data

        //- Minimum distance to patch
        scalar minPatchDistance_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        minPatchDistanceRefinement(const minPatchDistanceRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const minPatchDistanceRefinement&);


public:

    //- Runtime type information
    TypeName("minPatchDistanceRefinement");


    // Constructors

        //- Construct from components
        minPatchDistanceRefinement(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~minPatchDistanceRefinement();


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
