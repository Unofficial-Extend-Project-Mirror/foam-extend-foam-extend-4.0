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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace blockMatrixTools
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class BlockType>
void blockInsert
(
    const direction dir,
    const scalarField& x,
    Field<BlockType>& blockX
)
{
    forAll (x, i)
    {
        blockX[i](dir) = x[i];
    }
}


template<class BlockType>
void blockRetrieve
(
    const direction dir,
    scalarField& x,
    const Field<BlockType>& blockX
)
{
    forAll (x, i)
    {
        x[i] = blockX[i](dir);
    }
}


template<class BlockType>
void insertDiagSource
(
    const direction dir,
    const fvScalarMatrix& m,
    BlockLduMatrix<BlockType>& blockM,
    Field<BlockType>& blockB
)
{
    // Prepare the diagonal and source

    scalarField diag = m.diag();
    scalarField source = m.source();

    // Add boundary source contribution
    m.addBoundaryDiag(diag, 0);
    m.addBoundarySource(source, false);

    if (blockM.diag().activeType() == blockCoeffBase::UNALLOCATED)
    {
        blockM.diag().asScalar() = diag;
    }
    else if
    (
        blockM.diag().activeType() == blockCoeffBase::SCALAR
     || blockM.diag().activeType() == blockCoeffBase::LINEAR
    )
    {
        typename CoeffField<BlockType>::linearTypeField& blockDiag =
            blockM.diag().asLinear();

        forAll (diag, i)
        {
            blockDiag[i](dir) = diag[i];
        }
    }
    else if (blockM.diag().activeType() == blockCoeffBase::SQUARE)
    {
        typename CoeffField<BlockType>::squareTypeField& blockDiag =
            blockM.diag().asSquare();

        forAll (diag, i)
        {
            blockDiag[i](dir, dir) = diag[i];
        }
    }

    blockInsert(dir, source, blockB);
}


template<class BlockType>
void insertUpperLower
(
    const direction dir,
    const fvScalarMatrix& m,
    BlockLduMatrix<BlockType>& blockM
)
{
    if (m.diagonal())
    {
        // Matrix for insertion is diagonal-only: nothing to do
        return;
    }

    if (m.hasUpper())
    {
        const scalarField& upper = m.upper();

        if (blockM.upper().activeType() == blockCoeffBase::UNALLOCATED)
        {
            blockM.upper().asScalar() = upper;
        }
        else if
        (
            blockM.upper().activeType() == blockCoeffBase::SCALAR
         || blockM.upper().activeType() == blockCoeffBase::LINEAR
        )
        {
            typename CoeffField<BlockType>::linearTypeField& blockUpper =
                blockM.upper().asLinear();

            forAll (upper, i)
            {
                blockUpper[i](dir) = upper[i];
            }
        }
        else if (blockM.upper().activeType() == blockCoeffBase::SQUARE)
        {
            typename CoeffField<BlockType>::squareTypeField& blockUpper =
                blockM.upper().asSquare();

            forAll (upper, i)
            {
                blockUpper[i](dir, dir) = upper[i];
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void insertUpperLower\n"
            "(\n"
            "    const direction dir,\n"
            "    const fvScalarMatrix& m,\n"
            "    BlockLduMatrix<BlockType>& blockM\n"
            ")"
        )   << "Error in matrix insertion: problem with block structure"
            << abort(FatalError);
    }

    if (m.symmetric() && blockM.symmetric())
    {
        Info<< "Both m and blockM are symmetric: inserting only upper triangle"
            << endl;
    }
    else
    {
        // Either scalar or block matrix is asymmetric: insert lower triangle
        const scalarField& lower = m.lower();

        if (blockM.lower().activeType() == blockCoeffBase::UNALLOCATED)
        {
            blockM.lower().asScalar() = lower;
        }
        else if
        (
            blockM.lower().activeType() == blockCoeffBase::SCALAR
         || blockM.lower().activeType() == blockCoeffBase::LINEAR
        )
        {
            typename CoeffField<BlockType>::linearTypeField& blockLower =
                blockM.lower().asLinear();

            forAll (lower, i)
            {
                blockLower[i](dir) = lower[i];
            }
        }
        else if (blockM.lower().activeType() == blockCoeffBase::SQUARE)
        {
            typename CoeffField<BlockType>::squareTypeField& blockLower =
                blockM.lower().asSquare();

            forAll (lower, i)
            {
                blockLower[i](dir, dir) = lower[i];
            }
        }
    }
}


template<class BlockType>
void insertEquation
(
    const direction dir,
    const fvScalarMatrix& m,
    BlockLduMatrix<BlockType>& blockM,
    Field<BlockType>& blockX,
    Field<BlockType>& blockB
)
{
    insertDiagSource(dir, m, blockM, blockB);
    insertUpperLower(dir, m, blockM);
    blockInsert(dir, m.psi(), blockX);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockMatrixTools

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
