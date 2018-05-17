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
    Foam::minCellVolumeRefinement

Description
    Selection of refinement cells based on a minimum cell volume. All cells
    with volume larger than specified minimum get selected.

SourceFiles
    minCellVolumeRefinement.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef minCellVolumeRefinement_H
#define minCellVolumeRefinement_H

#include "refinementSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class minCellVolumeRefinement Declaration
\*---------------------------------------------------------------------------*/

class minCellVolumeRefinement
:
    public refinementSelection
{
    // Private data

        //- Minimum cell volume
        scalar minCellV_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        minCellVolumeRefinement(const minCellVolumeRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const minCellVolumeRefinement&);


public:

    //- Runtime type information
    TypeName("minCellVolumeRefinement");


    // Constructors

        //- Construct from components
        minCellVolumeRefinement(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~minCellVolumeRefinement();


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
