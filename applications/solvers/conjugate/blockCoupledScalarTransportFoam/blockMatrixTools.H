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

InNamespace
    Foam::

Description
    Block matrix insertion and retrieval tools

SourceFiles
    blockMatrixTools.C

\*---------------------------------------------------------------------------*/


#ifndef blockMatrixTools_H
#define blockMatrixTools_H

#include "blockLduMatrices.H"
#include "fvMatrices.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Namespace blockMatrixTools functions Declaration
\*---------------------------------------------------------------------------*/

namespace blockMatrixTools
{
    //- Insert field into block field
    template<class BlockType>
    void blockInsert
    (
        const direction dir,
        const scalarField& x,
        Field<BlockType>& blockX
    );

    //- Retrieve field from block field
    template<class BlockType>
    void blockRetrieve
    (
        const direction dir,
        scalarField& x,
        const Field<BlockType>& blockX
    );

    //- Insert matrix diagonal and source into the block system
    template<class BlockType>
    void insertDiagSource
    (
        const direction dir,
        const fvScalarMatrix& m,
        BlockLduMatrix<BlockType>& blockM,
        Field<BlockType>& blockB
    );

    // Insert matrix into the block system
    template<class Type, class BlockType>
    void insertUpperLower
    (
        const direction dir,
        const fvScalarMatrix& m,
        BlockLduMatrix<BlockType>& blockM
    );

    // Insert matrix into the block system
    template<class BlockType>
    void insertEquation
    (
        const direction dir,
        const fvScalarMatrix& m,
        BlockLduMatrix<BlockType>& blockM,
        Field<BlockType>& blockX,
        Field<BlockType>& blockB
    );

    // Update coupling of block system
    // Subtracts the block-coefficient coupling as specified by the user
    // from the source, leaving the implicit update given by linearisation
    template<class BlockType>
    void updateCoupling
    (
        BlockLduMatrix<BlockType>& blockM,
        Field<BlockType>& x,
        Field<BlockType>& b
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "blockMatrixTools.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
