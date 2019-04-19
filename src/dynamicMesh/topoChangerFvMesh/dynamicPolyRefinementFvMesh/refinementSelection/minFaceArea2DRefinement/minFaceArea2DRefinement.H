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
    Foam::minFaceArea2DRefinement

Description
    Selection of refinement cells based on a minimum face area of empty (or
    wedge) patches. All cells containing a face on empty (or wedge) patch larger
    than specified minimum area get selected.

    Useful as a 2D substitute for minCellVolume selection for wedge geometries.

SourceFiles
    minFaceArea2DRefinement.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef minFaceArea2DRefinement_H
#define minFaceArea2DRefinement_H

#include "refinementSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class minFaceArea2DRefinement Declaration
\*---------------------------------------------------------------------------*/

class minFaceArea2DRefinement
:
    public refinementSelection
{
    // Private data

        //- Minimum face area on empty (or wedge) patches
        scalar minFaceA_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        minFaceArea2DRefinement(const minFaceArea2DRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const minFaceArea2DRefinement&);


public:

    //- Runtime type information
    TypeName("minFaceArea2DRefinement");


    // Constructors

        //- Construct from components
        minFaceArea2DRefinement(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~minFaceArea2DRefinement();


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
