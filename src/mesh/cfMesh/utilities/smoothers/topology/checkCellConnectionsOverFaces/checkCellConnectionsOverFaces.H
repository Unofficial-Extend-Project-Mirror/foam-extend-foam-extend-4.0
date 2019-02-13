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
    checkCellConnectionsOverFaces

Description
    Checks cells if all cells can be visited over face neighbours. Deletes

SourceFiles
    checkCellConnectionsOverFaces.C

\*---------------------------------------------------------------------------*/

#ifndef checkCellConnectionsOverFaces_H
#define checkCellConnectionsOverFaces_H

#include "polyMeshGenModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class checkCellConnectionsOverFaces Declaration
\*---------------------------------------------------------------------------*/

class checkCellConnectionsOverFaces
{
    // Private data
        //- Reference to polyMeshGen
        polyMeshGen& mesh_;

        //- each cell is assigned a group
        //- cells which can be visited over face neighbours are assigned
        //- to the same group
        labelList cellGroup_;

        //- number of groups
        label nGroups_;

    // Private member functions
        //- decompose marked cells
        void findCellGroups();

        //- Disallow default bitwise copy construct
        checkCellConnectionsOverFaces(const checkCellConnectionsOverFaces&);

        //- Disallow default bitwise assignment
        void operator=(const checkCellConnectionsOverFaces&);

public:

    // Constructors

        //- Construct from polyMeshGen
        checkCellConnectionsOverFaces(polyMeshGen& mesh);

    // Destructor

        ~checkCellConnectionsOverFaces();

    // Member Functions
        //- check if all cells are connected as a single domain
        //- return true if the mesh has been changed
        bool checkCellGroups();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
