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
void blockAdd
(
    const direction dir,
    const scalarField& x,
    Field<BlockType>& blockX
)
{
    forAll (x, i)
    {
        blockX[i](dir) += x[i];
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

    // Insert block interface fields
    forAll(blockM.interfaces(), patchI)
    {
        if (blockM.interfaces().set(patchI))
        {
            // Couple upper and lower
            const scalarField& cUpper = m.boundaryCoeffs()[patchI];
            const scalarField& cLower = m.internalCoeffs()[patchI];

            if
            (
                blockM.coupleUpper()[patchI].activeType()
             == blockCoeffBase::UNALLOCATED
            )
            {
                blockM.coupleUpper()[patchI].asScalar() = cUpper;
                blockM.coupleLower()[patchI].asScalar() = cLower;
            }
            else if
            (
                blockM.coupleUpper()[patchI].activeType()
             == blockCoeffBase::SCALAR
             || blockM.coupleUpper()[patchI].activeType()
             == blockCoeffBase::LINEAR
            )
            {
                typename CoeffField<BlockType>::linearTypeField& blockUpper =
                    blockM.coupleUpper()[patchI].asLinear();

                typename CoeffField<BlockType>::linearTypeField& blockLower =
                    blockM.coupleLower()[patchI].asLinear();

                forAll (cUpper, i)
                {
                    blockUpper[i](dir) = cUpper[i];
                    blockLower[i](dir) = cLower[i];
                }
            }
            else if
            (
                blockM.coupleUpper()[patchI].activeType()
             == blockCoeffBase::SQUARE
            )
            {
                typename CoeffField<BlockType>::squareTypeField& blockUpper =
                    blockM.coupleUpper()[patchI].asSquare();

                typename CoeffField<BlockType>::squareTypeField& blockLower =
                    blockM.coupleLower()[patchI].asSquare();

                forAll (cUpper, i)
                {
                    blockUpper[i](dir, dir) = cUpper[i];
                    blockLower[i](dir, dir) = cLower[i];
                }
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


template<class BlockType>
void insertCouplingDiagSource
(
    const direction dirI,
    const direction dirJ,
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


    // Add off-diagonal block coefficients
    typename CoeffField<BlockType>::squareTypeField& blockDiag =
        blockM.diag().asSquare();

    // Set off-diagonal coefficient
    forAll (diag, i)
    {
        blockDiag[i](dirI, dirJ) = diag[i];
    }

    blockAdd(dirI, source, blockB);
}


template<class BlockType>
void insertCouplingUpperLower
(
    const direction dirI,
    const direction dirJ,
    const fvScalarMatrix& m,
    BlockLduMatrix<BlockType>& blockM
)
{
    if (m.diagonal())
    {
        // Matrix for insertion is diagonal-only: nothing to do
        return;
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

        typename CoeffField<BlockType>::squareTypeField& blockLower =
            blockM.lower().asSquare();

        forAll (lower, i)
        {
            blockLower[i](dirI, dirJ) = lower[i];
        }
    }

    if (m.hasUpper())
    {
        const scalarField& upper = m.upper();

        typename CoeffField<BlockType>::squareTypeField& blockUpper =
            blockM.upper().asSquare();

        forAll (upper, i)
        {
            blockUpper[i](dirI, dirJ) = upper[i];
        }
    }
    else
    {
        FatalErrorIn
        (
            "void insertCouplingUpperLower\n"
            "(\n"
            "    const direction dirI,\n"
            "    const direction dirJ,\n"
            "    const fvScalarMatrix& m,\n"
            "    BlockLduMatrix<BlockType>& blockM\n"
            ")"
        )   << "Error in matrix insertion: problem with block structure"
            << abort(FatalError);
    }

    // Insert block interface fields
    forAll(blockM.interfaces(), patchI)
    {
        if (blockM.interfaces().set(patchI))
        {
            // Couple upper and lower
            const scalarField& cUpper = m.boundaryCoeffs()[patchI];
            const scalarField& cLower = m.internalCoeffs()[patchI];

            typename CoeffField<BlockType>::squareTypeField& blockUpper =
                blockM.coupleUpper()[patchI].asSquare();

            typename CoeffField<BlockType>::squareTypeField& blockLower =
                blockM.coupleLower()[patchI].asSquare();

            forAll (cUpper, i)
            {
                blockUpper[i](dirI, dirJ) = cUpper[i];
                blockLower[i](dirI, dirJ) = cLower[i];
            }
        }
    }
}


template<class BlockType>
void insertCoupling
(
    const direction dirI,
    const direction dirJ,
    const fvScalarMatrix& m,
    BlockLduMatrix<BlockType>& blockM,
    Field<BlockType>& blockX,
    Field<BlockType>& blockB
)
{
    insertCouplingDiagSource(dirI, dirJ, m, blockM, blockB);
    insertCouplingUpperLower(dirI, dirJ, m, blockM);
}


template<class BlockType>
void updateSourceCoupling
(
    BlockLduMatrix<BlockType>& blockM,
    Field<BlockType>& x,
    Field<BlockType>& b
)
{
    // Eliminate off-diagonal block coefficients from the square diagonal
    // With this change, the segregated matrix can be assembled with complete
    // source terms and linearisation can be provided independently.
    // Once the linearisation coefficients are set (off-diagonal entries
    // in the square block matrix, they are multiplied by the current value
    // of the field and subtracted from the source term

    if (blockM.diag().activeType() == blockCoeffBase::SQUARE)
    {
        typename CoeffField<BlockType>::squareTypeField& blockDiag =
            blockM.diag().asSquare();

        typename CoeffField<BlockType>::linearTypeField lf(blockDiag.size());
        typename CoeffField<BlockType>::squareTypeField sf(blockDiag.size());

        // Expand and contract

        // Take out the diagonal entries from the square coefficient
        contractLinear(lf, blockDiag);

        // Expand the diagonal for full square, with zeroes in the off-diagonal
        expandLinear(sf, lf);

        // Subtract from the source the difference between the full block
        // diagonal and the diagonal terms only
        // Sign is the same as in the derivative
        b += (blockDiag - sf) & x;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockMatrixTools

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
