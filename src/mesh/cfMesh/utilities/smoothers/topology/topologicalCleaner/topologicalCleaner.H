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
    topologicalCleaner

Description
    Checks cells for invalid topologies. Decomposes bad cells into pyramids

SourceFiles
    topologicalCleaner.C

\*---------------------------------------------------------------------------*/

#ifndef topologicalCleaner_H
#define topologicalCleaner_H

#include "polyMeshGenModifier.H"
#include "boolList.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class topologicalCleaner Declaration
\*---------------------------------------------------------------------------*/

class topologicalCleaner
{
    // Private data
        polyMeshGen& mesh_;

        bool changed_;

        boolList decomposeCell_;

    // Private member functions

        //- check and fix incorrectly connected faces
        void checkNonConsecutiveBoundaryVertices();

        //- check and fix non-mappable faces
        void checkNonMappableFaces();

        //- check and fix non-mappable cells
        void checkNonMappableCells();

        //- decompose marked cells
        void decomposeCells();

        //- Disallow default bitwise copy construct
        topologicalCleaner(const topologicalCleaner&);

        //- Disallow default bitwise assignment
        void operator=(const topologicalCleaner&);

public:

    // Constructors

        //- Construct from polyMeshGen
        topologicalCleaner(polyMeshGen& mesh);


    // Destructor

        ~topologicalCleaner();


    // Member Functions
        //- clean topologically problematic cells
        bool cleanTopology();

        //- check and fix vertices where two or more groups of cells meet
        void checkInvalidConnectionsForVerticesCells
        (
            labelHashSet* irregularNodes = nullptr
        );

        //- check and fix vertices where two or more groups of faces meet
        void checkInvalidConnectionsForVerticesFaces
        (
            labelHashSet* irregularNodes = nullptr
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
