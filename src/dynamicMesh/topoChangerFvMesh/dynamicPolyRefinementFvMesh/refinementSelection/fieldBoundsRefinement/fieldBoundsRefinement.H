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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fieldBoundsRefinement Declaration
\*---------------------------------------------------------------------------*/

class fieldBoundsRefinement
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
        fieldBoundsRefinement
        (
            const dynamicPolyRefinementFvMesh& dynamicRefMesh,
            const polyhedralRefinement& pRef,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fieldBoundsRefinement();


    // Member Functions

        // Selection of refinement/unrefinement candidates

            //- Return transferable list of cells to refine
            const Xfer<labelList> refinementCellCandidates() const;

            //- Return transferable list of split points to unrefine
            const Xfer<labelList> unrefinementPointCandidates() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
