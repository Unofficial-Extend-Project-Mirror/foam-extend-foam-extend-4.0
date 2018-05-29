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
