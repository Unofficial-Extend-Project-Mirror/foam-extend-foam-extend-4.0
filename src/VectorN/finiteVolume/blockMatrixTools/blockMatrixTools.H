/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    // Field operations

        //- Insert field into block field
        template<class BlockType>
        void blockInsert
        (
            const direction dir,
            const scalarField& x,
            Field<BlockType>& blockX
        );

        //- Add field into block field
        template<class BlockType>
        void blockAdd
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


    // Diagonal block operations

        //- Insert matrix diagonal and source into the block system
        template<class BlockType>
        void insertDiagSource
        (
            const direction dir,
            const fvScalarMatrix& m,
            BlockLduMatrix<BlockType>& blockM,
            Field<BlockType>& blockB
        );

        // Insert upper and lower coefficients matrix into the block system
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


    // Coupling block operations

        //- Insert coupling matrix diagonal and source into the block system
        template<class BlockType>
        void insertCouplingDiagSource
        (
            const direction dirI,
            const direction dirJ,
            const fvScalarMatrix& m,
            BlockLduMatrix<BlockType>& blockM,
            Field<BlockType>& blockB
        );

        // Insert coupling matrix into the block system
        template<class Type, class BlockType>
        void insertCouplingUpperLower
        (
            const direction dirI,
            const direction dirJ,
            const fvScalarMatrix& m,
            BlockLduMatrix<BlockType>& blockM
        );

        // Insert coupling matrix into the block system
        template<class BlockType>
        void insertCoupling
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
        void updateSourceCoupling
        (
            BlockLduMatrix<BlockType>& blockM,
            Field<BlockType>& x,
            Field<BlockType>& b
        );

    //- Tools specific to pU coupled block matrix

        //- Insert field into block field to be solved for
        template<class blockType, class fieldType>
        void insertSolutionVector
        (
            const label loc,
            const Field<fieldType>& xSingle,
            Field<blockType>& xBlock
        );

        //- Retrieve solved field from block field
        template<class blockType, class fieldType>
        void retrieveSolution
        (
            const label loc,
            Field<fieldType>& xSingle,
            const Field<blockType>& xBlock
        );

        //- Insert matrix diagonal and source into the block system
        template<class blockType, class matrixType>
        void insertDiagSource
        (
            const label loc,
            fvMatrix<matrixType>& matrix,
            BlockLduMatrix<blockType>& A,
            Field<blockType>& b
        );

        // Insert upper and lower part into the block system
        template<class blockType, class matrixType>
        void insertUpperLower
        (
            const label loc,
            const fvMatrix<matrixType>& matrix,
            BlockLduMatrix<blockType>& A
        );

        // Update coupling coefficients in the block matrix
        template<class blockType, class matrixType>
        void updateCouplingCoeffs
        (
            const label loc,
            const fvMatrix<matrixType>& matrix,
            BlockLduMatrix<blockType>& A
        );

        // Insert matrix into the block system
        template<class blockType, class matrixType>
        void insertEquation
        (
            const label loc,
            fvMatrix<matrixType>& matrix,
            BlockLduMatrix<blockType>& A,
            Field<blockType>& x,
            Field<blockType>& b
        );

        // Insert block coupling
        template<class blockType1, class blockType2>
        void insertBlock
        (
            const label loc1,
            const label loc2,
            const BlockLduMatrix<blockType1>& blockMatrix,
            BlockLduMatrix<blockType2>& A,
            const bool incFirst
        );

        // Insert boundary contributions. These functions update the matrix coeffs
        // according to boundary conditions. Two functions are needed for
        // specialization if the field is scalar or vector. If the field is scalar,
        // it is assumed that the grad(p) was used to obtain the matrix, and the
        // boundary contributions to the source are vectors. If the field is
        // vector, assumed is the div(U) and the source is scalar.
        // Note parameter GeometricField<scalar...>, this is grad specialization.
        template<class blockType1, class blockType2>
        void insertBoundaryContributions
        (
            const label loc1,
            const label loc2,
            const BlockLduMatrix<blockType1>& blockMatrix,
            GeometricField<scalar, fvPatchField, volMesh>& psi,
            BlockLduMatrix<blockType2>& A,
            Field<blockType2>& b,
            const bool incFirst
        );

        // Note parameter GeometricField<vector...>, this is div specialization.
        template<class blockType1, class blockType2>
        void insertBoundaryContributions
        (
            const label loc1,
            const label loc2,
            const BlockLduMatrix<blockType1>& blockMatrix,
            GeometricField<vector, fvPatchField, volMesh>& psi,
            BlockLduMatrix<blockType2>& A,
            Field<blockType2>& b,
            const bool incFirst
        );

        // Insert existing block matrix (obtained by implicit grad/div operator)
        // into block system, updating the matrix coefficients from boundaries.
        // To be used ONLY for pressure velocity coupling terms, other scenarios
        // are not considered at the moment.
        template<class blockType1, class blockType2, class fieldType>
        void insertBlockCoupling
        (
            const label loc1,
            const label loc2,
            const BlockLduMatrix<blockType1>& blockMatrix,
            GeometricField<fieldType, fvPatchField, volMesh>& psi,
            BlockLduMatrix<blockType2>& A,
            Field<blockType2>& b,
            const bool incFirst
        );

        // Insert scalar equation coupling. Not tested: VV, 9/May/2014
        template<class blockType>
        void insertEquationCoupling
        (
            const label loc1,
            const label loc2,
            fvScalarMatrix& matrix,
            BlockLduMatrix<blockType>& A,
            Field<blockType> b
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
