/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    Foam::FECCellToFaceStencil

Description
    All cells connected via edge to face.

SourceFiles
    FECCellToFaceStencil.C

\*---------------------------------------------------------------------------*/

#ifndef FECCellToFaceStencil_H
#define FECCellToFaceStencil_H

#include "cellToFaceStencil.H"
#include "EdgeMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FECCellToFaceStencil Declaration
\*---------------------------------------------------------------------------*/

class FECCellToFaceStencil
:
    public cellToFaceStencil
{
    // Private Member Functions

        //- Calculates per edge the neighbour data (= edgeCells)
        void calcEdgeBoundaryData
        (
            const boolList& isValidBFace,
            const labelList& boundaryEdges,
            EdgeMap<labelList>& neiGlobal
        ) const;

        void calcFaceStencil(labelListList& faceStencil) const;

public:

    // Constructors

        //- Construct from all cells and boundary faces
        explicit FECCellToFaceStencil(const polyMesh&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
