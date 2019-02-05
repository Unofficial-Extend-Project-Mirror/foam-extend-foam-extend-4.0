/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    renameBoundaryPatches

Description
    Renames boundary patches based on user intput and sets correct patch types

SourceFiles
    renameBoundaryPatches.C

\*---------------------------------------------------------------------------*/

#ifndef renameBoundaryPatches_H
#define renameBoundaryPatches_H

#include "polyMeshGenModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class IOdictionary;

/*---------------------------------------------------------------------------*\
                           Class renameBoundaryPatches Declaration
\*---------------------------------------------------------------------------*/

class renameBoundaryPatches
{
    // Private data
        //- mesh
        polyMeshGen& mesh_;

        //- dictionary containing relevant information
        const IOdictionary& meshDict_;

    // Private member functions
        //- calculate new boundary and replace the existing one
        void calculateNewBoundary();

        //- find empty patches in 3D meshes and change them to wall
        void checkEmptyPatches();

        //- check if there exist any symmetry planes
        //- and make sure they get perfectly flat
        void checkSymmetryPlanes();

        //- Disallow default bitwise copy construct
        renameBoundaryPatches(const renameBoundaryPatches&);

        //- Disallow default bitwise assignment
        void operator=(const renameBoundaryPatches&);

public:

    // Constructors

        //- Construct from mesh and IOdictionary
        renameBoundaryPatches
        (
            polyMeshGen& mesh,
            const IOdictionary& meshDict,
            const bool allowEmptyPatches = false
        );

    // Destructor

        ~renameBoundaryPatches();

    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
