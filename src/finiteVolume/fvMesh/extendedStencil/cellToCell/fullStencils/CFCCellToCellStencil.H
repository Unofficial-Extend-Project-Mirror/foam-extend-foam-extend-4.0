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
    Foam::CFCCellToCellStencil

Description

SourceFiles
    CFCCellToCellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef CFCCellToCellStencil_H
#define CFCCellToCellStencil_H

#include "cellToCellStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class CFCCellToCellStencil Declaration
\*---------------------------------------------------------------------------*/

class CFCCellToCellStencil
:
    public cellToCellStencil
{
    // Private Member Functions

        void calcFaceBoundaryData(labelList& neiGlobal) const;

        void calcCellStencil(labelListList& globalCellCells) const;

        //- Disallow default bitwise copy construct
        CFCCellToCellStencil(const CFCCellToCellStencil&);

        //- Disallow default bitwise assignment
        void operator=(const CFCCellToCellStencil&);

public:

    // Constructors

        //- Construct from mesh
        explicit CFCCellToCellStencil(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
