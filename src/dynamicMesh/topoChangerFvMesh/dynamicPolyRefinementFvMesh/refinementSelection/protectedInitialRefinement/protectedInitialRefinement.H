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
    Foam::protectedInitialRefinement

Description
    This refinement strategy is just a protection to avoid unrefining cells that
    are refined during the snappyHexMesh meshing process. It's intended use is
    alongside some other refinement strategy as a part of a composite scheme.
    Therefore: refinementCellCandidates returns all cells as refinement
    candidates, while unrefinementPointCandidates returns all points that don't
    touch any of the initially refined cells.

    Since snappyHexMesh writes down cellLevel and pointLevel, cells refined in
    during meshing could be unrefined, but this is not what we would like to
    achieve (after all, if we have included refinement regions in meshing, then
    we probably had a good reason for it).

SourceFiles
    protectedInitialRefinement.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef protectedInitialRefinement_H
#define protectedInitialRefinement_H

#include "refinementSelection.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class protectedInitialRefinement Declaration
\*---------------------------------------------------------------------------*/

class protectedInitialRefinement
:
    public refinementSelection
{
    // Private data

        //- Initial cellLevel as volScalarField in order to simplify data
        //  transfer for dynamic load balancing
        volScalarField initialCellLevel_;

        //- Minimum refinement level to protect
        scalar minProtectedLevel_;


    // Private Member Functions

        // Copy control

            //- Disallow default bitwise copy construct
            protectedInitialRefinement(const protectedInitialRefinement&);

            //- Disallow default bitwise assignment
            void operator=(const protectedInitialRefinement&);


        // Constructor helper

            //- Create and return volScalarField from cellLevel
            tmp<volScalarField> cellLevelAsField() const;


public:

    //- Runtime type information
    TypeName("protectedInitialRefinement");


    // Constructors

        //- Construct from components
        protectedInitialRefinement(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~protectedInitialRefinement();


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
