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
    Foam::refinementSelection

Description
    An abstract base class for providing an interface for refinement selection
    algorithms. The interface provides two functionalities:

    1. Selects candidate cells to refine based on chosen criteria through
       refinementCellCandidates() member function,
    2. Selects candidate split points to unrefine based on chosen criteria
       through unrefinementPointCandidates() member function.

    Note: Here, we do not check that the unrefinement point candidates do not
    clash the points of cells marked as candidates for refinement. This is taken
    care of in polyhedralRefinement class.

SourceFiles
    refinementSelection.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef refinementSelection_H
#define refinementSelection_H

#include "dictionary.H"
#include "Xfer.H"
#include "labelList.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dynamicPolyRefinementFvMesh;
class polyhedralRefinement;

/*---------------------------------------------------------------------------*\
                   Class refinementSelection Declaration
\*---------------------------------------------------------------------------*/

class refinementSelection
{
    // Private data

        //- Const reference to dynamicPolyRefinementFvMesh
        const dynamicPolyRefinementFvMesh& dynamicRefMesh_;

        //- Const reference to polyhedralRefinement mesh modifier
        const polyhedralRefinement& pRef_;

        //- Refinement selection dictionary (subdictionary of
        //  dynamicPolyRefinementFvMeshCoeffs dictionary)
        const dictionary coeffDict_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        refinementSelection(const refinementSelection&);

        //- Disallow default bitwise assignment
        void operator=(const refinementSelection&);


protected:

    // Protected member functions

        // Access functions for derived classes

            //- Const access to the dynamic refinement mesh
            const dynamicPolyRefinementFvMesh& dynamicRefMesh() const
            {
                return dynamicRefMesh_;
            }

            //- Const access to the polyhedral refinement mesh modifier
            const polyhedralRefinement& polyRef() const
            {
                return polyRef_;
            }

            //- Const access to the coefficient dictionary
            const dictionary& coeffDict() const
            {
                return coeffDict_;
            }


public:

    //- Runtime type information
    TypeName("refinementSelection");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            refinementSelection,
            dictionary,
            (
                const dynamicPolyRefinementFvMesh& dynamicRefMesh,
                const polyhedralRefinement& pRef,
                const dictionary& dict
            ),
            (dynamicRefMesh, pRef, dict)
        );


    // Constructors

        //- Construct from components
        refinementSelection
        (
            const dynamicPolyRefinementFvMesh& dynamicRefMesh,
            const polyhedralRefinement& pRef,
            const dictionary& dict
        );


    // Selectors

        //- Return an autoPtr to the selected refinementSelection
        static autoPtr<refinementSelection> New
        (
            const dynamicPolyRefinementFvMesh& dynamicRefMesh,
            const polyhedralRefinement& pRef,
            const dictionary& dict
        );


    //- Destructor
    virtual ~refinementSelection();


    // Member Functions

        // Selection of refinement/unrefinement candidates

            //- Return transferable list of refinement cell candidates
            const Xfer<labelList> refinementCellCandidates() const = 0;

            //- Return transferable list of unrefinement split point candidates
            const Xfer<labelList> unrefinementPointCandidates() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
