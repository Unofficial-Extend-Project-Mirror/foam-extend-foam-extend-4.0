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
    Foam::fieldBoundsRefinement

Description
    Selection of refinement cells based on a given scalar field bounds:
    if a given (scalar) field is between lowerBound and upperBound, cell will
    be marked for refinement.

SourceFiles
    fieldBoundsRefinement.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef fieldBoundsRefinement_H
#define fieldBoundsRefinement_H

#include "refinementSelection.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fieldBoundsRefinement Declaration
\*---------------------------------------------------------------------------*/

class fieldBoundsRefinement
:
    public refinementSelection
{
    // Private data

        //- Name of the field to refine
        word fieldName_;

        //- Lower bound for the field
        scalar lowerBound_;

        //- Upper bound for the field
        scalar upperBound_;

        //- Whether to use cell-point-cell interpolation to smooth out the field
        //  before selection. Switched off by default
        Switch cellPointCellSmoothing_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fieldBoundsRefinement(const fieldBoundsRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const fieldBoundsRefinement&);


public:

    //- Runtime type information
    TypeName("fieldBoundsRefinement");


    // Constructors

        //- Construct from components
        fieldBoundsRefinement(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~fieldBoundsRefinement();


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
