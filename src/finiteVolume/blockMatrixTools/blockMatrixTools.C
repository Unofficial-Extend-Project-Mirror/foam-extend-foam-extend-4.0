/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

template<class blockType, class fieldType>
void insertSolutionVector
(
    const label loc,
    const Field<fieldType>& xSingle,
    Field<blockType>& xBlock
)
{
    // Get number of field components and local copy of location, for
    // consistency with member functions where locations need to be reset.
    const label nCmpts = pTraits<fieldType>::nComponents;
    label localLoc = loc;

    for (label cmptI = 0; cmptI < nCmpts; cmptI++)
    {
        scalarField xSingleCurr(xSingle.component(cmptI));

        forAll (xSingleCurr, cellI)
        {
            xBlock[cellI](localLoc) = xSingleCurr[cellI];
        }

        localLoc++;
    }
}


template<class blockType, class fieldType>
void retrieveSolution
(
    const label loc,
    Field<fieldType>& xSingle,
    const Field<blockType>& xBlock
)
{
    const label nCmpts = pTraits<fieldType>::nComponents;
    label localLoc = loc;

    for (label cmptI = 0; cmptI < nCmpts; cmptI++)
    {
        scalarField xSingleCurr(xSingle.component(cmptI));

        forAll (xSingleCurr, cellI)
        {
            xSingleCurr[cellI] = xBlock[cellI](localLoc);
        }

        xSingle.replace(cmptI, xSingleCurr);

        localLoc++;
    }
}


template<class blockType, class matrixType>
void insertDiagSource
(
    const label loc,
    fvMatrix<matrixType>& matrix,
    BlockLduMatrix<blockType>& A,
    Field<blockType>& b
)
{
    matrix.completeAssembly();

    // Save a copy for different components
    scalarField& diag = matrix.diag();
    scalarField saveDiag(diag);

    // Add source boundary contribution
    Field<matrixType>& source = matrix.source();
    matrix.addBoundarySource(source, false);

    const label nCmpts = pTraits<matrixType>::nComponents;
    label localLoc = loc;

    if
    (
        // This is needed if the matrixType is <vector>, then you need to grab
        // coeffs as linear. Consider doing a matrixType check also.
        // VV, 17/March/2014
        A.diag().activeType() != blockCoeffBase::SQUARE
    )
    {
        typename CoeffField<blockType>::linearTypeField& blockDiag =
            A.diag().asLinear();

        for (label cmptI = 0; cmptI < nCmpts; cmptI++)
        {
            matrix.addBoundaryDiag(diag, cmptI);
            scalarField sourceCmpt(source.component(cmptI));

//            FieldField<Field, scalar> bouCoeffsCmpt
//            (
//                matrix.boundaryCoeffs().component(cmptI)
//            );

            // Possible problem for coupled non-aligned boundaries.
            // VV, 14/May/2014.
//            matrix.correctImplicitBoundarySource
//            (
//                bouCoeffsCmpt,
//                sourceCmpt,
//                cmptI
//            );

            forAll (diag, cellI)
            {
                blockDiag[cellI](localLoc) = diag[cellI];
                b[cellI](localLoc) += sourceCmpt[cellI];
            }

            localLoc++;

            // Reset diagonal
            diag = saveDiag;
        }
    }
    else if (A.diag().activeType() == blockCoeffBase::SQUARE)
    {
        typename CoeffField<blockType>::squareTypeField& blockDiag =
            A.diag().asSquare();

        for (label cmptI = 0; cmptI < nCmpts; cmptI++)
        {
            matrix.addBoundaryDiag(diag, cmptI);
            scalarField sourceCmpt(source.component(cmptI));

//            FieldField<Field, scalar> bouCoeffsCmpt
//            (
//                matrix.boundaryCoeffs().component(cmptI)
//            );

//            matrix.correctImplicitBoundarySource
//            (
//                bouCoeffsCmpt,
//                sourceCmpt,
//                cmptI
//            );

            forAll (diag, cellI)
            {
                blockDiag[cellI](localLoc, localLoc) = diag[cellI];
                b[cellI](localLoc) += sourceCmpt[cellI];
            }

            localLoc++;

            // Reset diagonal
            diag = saveDiag;
        }
    }
}


template<class blockType, class matrixType>
void insertUpperLower
(
    const label loc,
    const fvMatrix<matrixType>& matrix,
    BlockLduMatrix<blockType>& A
)
{
    if (matrix.diagonal())
    {
        // Matrix for insertion is diagonal-only: nothing to do
        return;
    }

    const label nCmpts = pTraits<matrixType>::nComponents;
    label localLoc = loc;

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();

        if (A.upper().activeType() == blockCoeffBase::UNALLOCATED)
        {
            A.upper().asScalar() = upper;
        }
        else if
        (
            A.upper().activeType() == blockCoeffBase::SCALAR
         || A.upper().activeType() == blockCoeffBase::LINEAR
        )
        {
            typename CoeffField<blockType>::linearTypeField& blockUpper =
                A.upper().asLinear();

            for (label cmptI = 0; cmptI < nCmpts; cmptI++)
            {
                forAll (upper, faceI)
                {
                    blockUpper[faceI](localLoc) = upper[faceI];
                }

                localLoc++;
            }
        }
        else if (A.upper().activeType() == blockCoeffBase::SQUARE)
        {
            typename CoeffField<blockType>::squareTypeField& blockUpper =
                A.upper().asSquare();

            for (label cmptI = 0; cmptI < nCmpts; cmptI++)
            {
                forAll (upper, faceI)
                {
                    blockUpper[faceI](localLoc, localLoc) = upper[faceI];
                }

                localLoc++;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void insertUpperLower\n"
            "(\n"
            "    const location loc,\n"
            "    const fvMatrix<matrixType>& matrix,\n"
            "    BlockLduMatrix<blockType>& A\n"
            ")"
        )   << "Error in matrix insertion: problem with block structure."
            << abort(FatalError);
    }

    if (matrix.symmetric() && A.symmetric())
    {
        Info<< "Both matrices are symmetric: inserting only upper triangle"
            << endl;
    }
    else
    {
        // Reset localLoc
        localLoc = loc;

        // Either scalar or block matrix is asymmetric: insert lower triangle
        const scalarField& lower = matrix.lower();

        if (A.lower().activeType() == blockCoeffBase::UNALLOCATED)
        {
            A.lower().asScalar() = lower;
        }
        else if
        (
            A.lower().activeType() == blockCoeffBase::SCALAR
         || A.lower().activeType() == blockCoeffBase::LINEAR
        )
        {
            typename CoeffField<blockType>::linearTypeField& blockLower =
                A.lower().asLinear();

            for (label cmptI = 0; cmptI < nCmpts; cmptI++)
            {
                forAll (lower, faceI)
                {
                    blockLower[faceI](localLoc) = lower[faceI];
                }

                localLoc++;
            }
        }
        else if (A.lower().activeType() == blockCoeffBase::SQUARE)
        {
            typename CoeffField<blockType>::squareTypeField& blockLower =
                A.lower().asSquare();

            for (label cmptI = 0; cmptI < nCmpts; cmptI++)
            {
                forAll (lower, faceI)
                {
                    blockLower[faceI](localLoc, localLoc) = lower[faceI];
                }

                localLoc++;
            }
        }
    }
}


template<class blockType, class matrixType>
void updateCouplingCoeffs
(
    const label loc,
    const fvMatrix<matrixType>& matrix,
    BlockLduMatrix<blockType>& A
)
{
    const label nCmpts = pTraits<matrixType>::nComponents;
    label localLoc = loc;

    const GeometricField<matrixType, fvPatchField, volMesh>& psi = matrix.psi();
    forAll(psi.boundaryField(), patchI)
    {
        const fvPatchField<matrixType>& pf = psi.boundaryField()[patchI];
        const fvPatch& patch = pf.patch();

        if (patch.coupled())
        {
            const Field<matrixType>& icp = matrix.internalCoeffs()[patchI];
            const Field<matrixType>& bcp = matrix.boundaryCoeffs()[patchI];

            if (A.coupleUpper()[patchI].activeType() != blockCoeffBase::SQUARE)
            {
                typename CoeffField<blockType>::linearTypeField& pcoupleUpper =
                    A.coupleUpper()[patchI].asLinear();
                typename CoeffField<blockType>::linearTypeField& pcoupleLower =
                    A.coupleLower()[patchI].asLinear();

                for (label cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    scalarField icpCmpt = icp.component(cmptI);
                    scalarField bcpCmpt = bcp.component(cmptI);

                    forAll(pf, faceI)
                    {
                        pcoupleUpper[faceI](localLoc) = bcpCmpt[faceI];
                        pcoupleLower[faceI](localLoc) = icpCmpt[faceI];
                    }

                    localLoc++;
                }

                localLoc = loc;
            }
            else if
            (
                A.coupleUpper()[patchI].activeType() == blockCoeffBase::SQUARE
            )
            {
                typename CoeffField<blockType>::squareTypeField& pcoupleUpper =
                    A.coupleUpper()[patchI].asSquare();
                typename CoeffField<blockType>::squareTypeField& pcoupleLower =
                    A.coupleLower()[patchI].asSquare();

                for (label cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    scalarField icpCmpt = icp.component(cmptI);
                    scalarField bcpCmpt = bcp.component(cmptI);

                    forAll(pf, faceI)
                    {
                        pcoupleUpper[faceI](localLoc, localLoc) =
                            bcpCmpt[faceI];
                        pcoupleLower[faceI](localLoc, localLoc) =
                            icpCmpt[faceI];
                    }

                    localLoc++;
                }

                localLoc = loc;
            }
        }
    }
}


template<class blockType, class matrixType>
void insertEquation
(
    const label loc,
    fvMatrix<matrixType>& matrix,
    BlockLduMatrix<blockType>& A,
    Field<blockType>& x,
    Field<blockType>& b
)
{
    insertDiagSource(loc, matrix, A, b);
    insertUpperLower(loc, matrix, A);
    updateCouplingCoeffs(loc, matrix, A);
    insertSolutionVector(loc, matrix.psi().internalField(), x);
}


template<class blockType1, class fieldType, class blockType2>
void insertBlock
(
    const label loc1,
    const label loc2,
    const BlockLduSystem<blockType1, fieldType>& blockSystem,
    BlockLduMatrix<blockType2>& A,
    const bool incFirst
)
{
    // Sanity checks
    {
        const label blockMatrixSize = pTraits<blockType1>::nComponents;
        const label blockMatrixASize = pTraits<blockType2>::nComponents;

        if (blockMatrixSize > blockMatrixASize)
        {
            FatalErrorIn
            (
                "void insertBlock\n"
                "(\n"
                "    const label loc1,\n"
                "    const label loc2,\n"
                "    BlockLduMatrix<blockType1>& blockMatrix,\n"
                "    BlockLduMatrix<blockType2>& A\n"
                ")"
            )   << "Trying to insert a block matrix into smaller one."
                << abort(FatalError);
        }

        if (loc1 == loc2)
        {
            FatalErrorIn
            (
                "void insertCoupling\n"
                "(\n"
                "    const label loc1,\n"
                "    const label loc2,\n"
                "    BlockLduMatrix<blockType1>& blockMatrix,\n"
                "    BlockLduMatrix<blockType2>& A\n"
                ")"
            )   << "Trying to insert coupling in the position where equation "
                << "should be, since loc1 = loc2. Try using insertEquatiion "
                << "member function."
                << abort(FatalError);
        }
    }

    const label nCmpts = pTraits<blockType1>::nComponents;
    label localLoc1 = loc1;
    label localLoc2 = loc2;

    // Get references to ldu fields of blockMatrix always as linear
    const typename CoeffField<blockType1>::linearTypeField& bmd =
        blockSystem.diag().asLinear();
    const typename CoeffField<blockType1>::linearTypeField& bmu =
        blockSystem.upper().asLinear();
    const typename CoeffField<blockType1>::linearTypeField& bml =
        blockSystem.lower().asLinear();

    // Get references to ldu fields of A matrix always as square
    typename CoeffField<blockType2>::squareTypeField& blockDiag =
         A.diag().asSquare();
    typename CoeffField<blockType2>::squareTypeField& blockUpper =
         A.upper().asSquare();
    typename CoeffField<blockType2>::squareTypeField& blockLower =
         A.lower().asSquare();

    // Insert blockMatrix that represents coupling into larger system matrix
    for (label cmptI = 0; cmptI < nCmpts; cmptI++)
    {
        forAll(bmd, cellI)
        {
            blockDiag[cellI](localLoc1, localLoc2) +=
                bmd[cellI].component(cmptI);
        }

        forAll(bmu, faceI)
        {
            blockUpper[faceI](localLoc1, localLoc2) +=
                bmu[faceI].component(cmptI);
            blockLower[faceI](localLoc1, localLoc2) +=
                bml[faceI].component(cmptI);
        }

        if (incFirst)
        {
            localLoc1++;
        }
        else
        {
            localLoc2++;
        }
    }
}


template<class blockType1, class blockType2, class fieldType>
void insertBoundaryContributions
(
    const label loc1,
    const label loc2,
    const BlockLduSystem<blockType1, fieldType>& blockSystem,
    BlockLduMatrix<blockType2>& A,
    Field<blockType2>& b,
    const bool incFirst
)
{
    // Need to get reference to fvMesh instead of lduMesh
    const fvBoundaryMesh& bmesh = refCast<const fvMesh>(A.mesh()).boundary();

    const label nCmpts = pTraits<blockType1>::nComponents;
    const label nSrcCmpts = pTraits<fieldType>::nComponents;
    label localLoc1 = loc1;
    label localLoc2 = loc2;

    const Field<fieldType>& source = blockSystem.source();

    // Insert source from block system to rhs
    for (label cmptI = 0; cmptI < nSrcCmpts; cmptI++)
    {
        scalarField sourceCmpt(source.component(cmptI));

        forAll(b, cellI)
        {
            b[cellI](localLoc1) += sourceCmpt[cellI];
        }

        if (incFirst)
        {
            localLoc1++;
        }
        else
        {
            localLoc2++;
        }
    }

    // Reset local locations for coupling contributions
    localLoc1 = loc1;
    localLoc2 = loc2;

    // Insert coupling contributions into block matrix
    forAll(bmesh, patchI)
    {
        if (bmesh[patchI].coupled())
        {
            typename CoeffField<blockType2>::squareTypeField& pcoupleUpper =
                A.coupleUpper()[patchI].asSquare();
            typename CoeffField<blockType2>::squareTypeField& pcoupleLower =
                A.coupleLower()[patchI].asSquare();

            const typename CoeffField<blockType1>::linearTypeField& bmcu =
                blockSystem.coupleUpper()[patchI].asLinear();
            const typename CoeffField<blockType1>::linearTypeField& bmcl =
                blockSystem.coupleLower()[patchI].asLinear();

                for (label cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    forAll(bmcu, faceI)
                    {
                        pcoupleUpper[faceI](localLoc1, localLoc2) +=
                              bmcu[faceI].component(cmptI);
                        pcoupleLower[faceI](localLoc1, localLoc2) +=
                              bmcl[faceI].component(cmptI);
                    }

                    if (incFirst)
                    {
                        localLoc1++;
                    }
                    else
                    {
                        localLoc2++;
                    }
                }

            // Reset local locations for other patches
            localLoc1 = loc1;
            localLoc2 = loc2;
        }
    }
}


template<class blockType1, class blockType2, class fieldType>
void insertBlockCoupling
(
    const label loc1,
    const label loc2,
    const BlockLduSystem<blockType1, fieldType>& blockSystem,
    BlockLduMatrix<blockType2>& A,
    Field<blockType2>& b,
    const bool incFirst
)
{
    insertBlock(loc1, loc2, blockSystem, A, incFirst);
    insertBoundaryContributions(loc1, loc2, blockSystem, A, b, incFirst);
}


template<class blockType>
void insertEquationCoupling
(
    const label loc1,
    const label loc2,
    fvScalarMatrix& matrix,
    BlockLduMatrix<blockType>& A,
    Field<blockType>& b
)
{
    // Sanity check
    if (loc1 == loc2)
    {
        FatalErrorIn
        (
            "void insertEquationCoupling\n"
            "(\n"
            "    const label loc1,\n"
            "    const label loc2,\n"
            "    fvScalarMatrix& matrix\n"
            "    BlockLduMatrix<blockType>& A\n"
            "    Field<blockType>& b\n"
            ")"
        )   << "Trying to insert coupling in the position where equation "
            << "should be since loc1 = loc2. Try using insertEquatiion "
            << "member function."
            << abort(FatalError);
    }

    // Get references to fvScalarMatrix fields, updating boundary contributions
    scalarField& diag = matrix.D();
    scalarField& source = matrix.source();
    matrix.addBoundarySource(source, false);

    const scalarField& upper = matrix.upper();
    const scalarField& lower = matrix.lower();

    // Get references to ldu fields of A matrix always as square
    typename CoeffField<blockType>::squareTypeField& blockDiag =
         A.diag().asSquare();
    typename CoeffField<blockType>::squareTypeField& blockUpper =
         A.upper().asSquare();
    typename CoeffField<blockType>::squareTypeField& blockLower =
         A.lower().asSquare();

    forAll(diag, cellI)
    {
        blockDiag[cellI](loc1, loc2) += diag[cellI];
        b[cellI](loc1) += source[cellI];
    }

    forAll(upper, faceI)
    {
        blockUpper[faceI](loc1, loc2) += upper[faceI];
        blockLower[faceI](loc1, loc2) += lower[faceI];
    }

    // Update coupling contributions
    const volScalarField& psi = matrix.psi();
    forAll(psi.boundaryField(), patchI)
    {
        const fvPatchScalarField& pf = psi.boundaryField()[patchI];
        const fvPatch& patch = pf.patch();

        if (patch.coupled())
        {
            const scalarField& icp = matrix.internalCoeffs()[patchI];
            const scalarField& bcp = matrix.boundaryCoeffs()[patchI];

            typename CoeffField<blockType>::squareTypeField& pcoupleUpper =
                A.coupleUpper()[patchI].asSquare();
            typename CoeffField<blockType>::squareTypeField& pcoupleLower =
                A.coupleLower()[patchI].asSquare();

             forAll(pf, faceI)
             {
                pcoupleUpper[faceI](loc1, loc2) = bcp[faceI];
                pcoupleLower[faceI](loc1, loc2) = icp[faceI];
             }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockMatrixTools

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
