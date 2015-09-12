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

#include "fvBlockMatrix.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class fieldType>
void fvBlockMatrix<Type>::insertSolutionVector
(
    const direction dir,
    const Field<fieldType>& xSingle
)
{
    // Get number of field components and local copy of direction, for
    // consistency with member functions where directions need to be reset.
    const direction nCmpts = pTraits<fieldType>::nComponents;
    direction localDir = dir;

    Field<Type>& psiIn = psi_.internalField();

    for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
    {
        scalarField xSingleCurr(xSingle.component(cmptI));

        forAll (xSingleCurr, cellI)
        {
            psiIn[cellI](localDir) = xSingleCurr[cellI];
        }

        localDir++;
    }
}


template<class Type>
template<class matrixType>
void fvBlockMatrix<Type>::insertDiagSource
(
    const direction dir,
    fvMatrix<matrixType>& matrix
)
{
    matrix.completeAssembly();

    // Save a copy for different components
    scalarField& diag = matrix.diag();
    scalarField saveDiag(diag);

    // Add source boundary contribution
    Field<matrixType>& source = matrix.source();
    matrix.addBoundarySource(source, false);

    const direction nCmpts = pTraits<matrixType>::nComponents;
    direction localDir = dir;

    // Get reference to this source field of block system
    Field<Type>& b = this->source();

    if
    (
        // This is needed if the matrixType is <vector>, then you need to grab
        // coeffs as linear. Consider doing a matrixType check also.
        // VV, 17/March/2014
        this->diag().activeType() != blockCoeffBase::SQUARE
    )
    {
        typename CoeffField<Type>::linearTypeField& blockDiag =
            this->diag().asLinear();

        for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
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
                blockDiag[cellI](localDir) = diag[cellI];
                b[cellI](localDir) += sourceCmpt[cellI];
            }

            localDir++;

            // Reset diagonal
            diag = saveDiag;
        }
    }
    else if (this->diag().activeType() == blockCoeffBase::SQUARE)
    {
        typename CoeffField<Type>::squareTypeField& blockDiag =
            this->diag().asSquare();

        for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
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
                blockDiag[cellI](localDir, localDir) = diag[cellI];
                b[cellI](localDir) += sourceCmpt[cellI];
            }

            localDir++;

            // Reset diagonal
            diag = saveDiag;
        }
    }
}


template<class Type>
template<class matrixType>
void fvBlockMatrix<Type>::insertUpperLower
(
    const direction dir,
    const fvMatrix<matrixType>& matrix
)
{
    if (matrix.diagonal())
    {
        // Matrix for insertion is diagonal-only: nothing to do
        return;
    }

    const direction nCmpts = pTraits<matrixType>::nComponents;
    direction localDir = dir;

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();

        if (this->upper().activeType() == blockCoeffBase::UNALLOCATED)
        {
            this->upper().asScalar() = upper;
        }
        else if
        (
            this->upper().activeType() == blockCoeffBase::SCALAR
         || this->upper().activeType() == blockCoeffBase::LINEAR
        )
        {
            typename CoeffField<Type>::linearTypeField& blockUpper =
                this->upper().asLinear();

            for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
            {
                forAll (upper, faceI)
                {
                    blockUpper[faceI](localDir) = upper[faceI];
                }

                localDir++;
            }
        }
        else if (this->upper().activeType() == blockCoeffBase::SQUARE)
        {
            typename CoeffField<Type>::squareTypeField& blockUpper =
                this->upper().asSquare();

            for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
            {
                forAll (upper, faceI)
                {
                    blockUpper[faceI](localDir, localDir) = upper[faceI];
                }

                localDir++;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void fvBlockMatrix<Type, matrixType>::insertUpperLower\n"
            "(\n"
            "    const direction dir,\n"
            "    const fvMatrix<matrixType>& matrix\n"
            ")"
        )   << "Error in matrix insertion: problem with block structure."
            << abort(FatalError);
    }

    if (matrix.symmetric() && this->symmetric())
    {
        Info<< "Both matrices are symmetric: inserting only upper triangle"
            << endl;
    }
    else
    {
        // Reset localDir
        localDir = dir;

        // Either scalar or block matrix is asymmetric: insert lower triangle
        const scalarField& lower = matrix.lower();

        if (this->lower().activeType() == blockCoeffBase::UNALLOCATED)
        {
            this->lower().asScalar() = lower;
        }
        else if
        (
            this->lower().activeType() == blockCoeffBase::SCALAR
         || this->lower().activeType() == blockCoeffBase::LINEAR
        )
        {
            typename CoeffField<Type>::linearTypeField& blockLower =
                this->lower().asLinear();

            for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
            {
                forAll (lower, faceI)
                {
                    blockLower[faceI](localDir) = lower[faceI];
                }

                localDir++;
            }
        }
        else if (this->lower().activeType() == blockCoeffBase::SQUARE)
        {
            typename CoeffField<Type>::squareTypeField& blockLower =
                this->lower().asSquare();

            for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
            {
                forAll (lower, faceI)
                {
                    blockLower[faceI](localDir, localDir) = lower[faceI];
                }

                localDir++;
            }
        }
    }
}


template<class Type>
template<class matrixType>
void fvBlockMatrix<Type>::updateCouplingCoeffs
(
    const direction dir,
    const fvMatrix<matrixType>& matrix
)
{
    const direction nCmpts = pTraits<matrixType>::nComponents;
    direction localDir = dir;

    const GeometricField<matrixType, fvPatchField, volMesh>& psi =
        matrix.psi();

    forAll(psi.boundaryField(), patchI)
    {
        const fvPatchField<matrixType>& pf = psi.boundaryField()[patchI];
        const fvPatch& patch = pf.patch();

        if (patch.coupled())
        {
            const Field<matrixType>& icp = matrix.internalCoeffs()[patchI];
            const Field<matrixType>& bcp = matrix.boundaryCoeffs()[patchI];

            if
            (
                this->coupleUpper()[patchI].activeType()
             != blockCoeffBase::SQUARE
            )
            {
                typename CoeffField<Type>::linearTypeField& pcoupleUpper =
                    this->coupleUpper()[patchI].asLinear();
                typename CoeffField<Type>::linearTypeField& pcoupleLower =
                    this->coupleLower()[patchI].asLinear();

                for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    scalarField icpCmpt = icp.component(cmptI);
                    scalarField bcpCmpt = bcp.component(cmptI);

                    forAll(pf, faceI)
                    {
                        pcoupleUpper[faceI](localDir) = bcpCmpt[faceI];
                        pcoupleLower[faceI](localDir) = icpCmpt[faceI];
                    }

                    localDir++;
                }

                localDir = dir;
            }
            else if
            (
                this->coupleUpper()[patchI].activeType()
             == blockCoeffBase::SQUARE
            )
            {
                typename CoeffField<Type>::squareTypeField& pcoupleUpper =
                    this->coupleUpper()[patchI].asSquare();
                typename CoeffField<Type>::squareTypeField& pcoupleLower =
                    this->coupleLower()[patchI].asSquare();

                for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    scalarField icpCmpt = icp.component(cmptI);
                    scalarField bcpCmpt = bcp.component(cmptI);

                    forAll(pf, faceI)
                    {
                        pcoupleUpper[faceI](localDir, localDir) =
                            bcpCmpt[faceI];
                        pcoupleLower[faceI](localDir, localDir) =
                            icpCmpt[faceI];
                    }

                    localDir++;
                }

                localDir = dir;
            }
        }
    }
}


template<class Type>
template<class blockType, class fieldType>
void Foam::fvBlockMatrix<Type>::insertBlock
(
    const direction dirI,
    const direction dirJ,
    const BlockLduSystem<blockType, fieldType>& blockSystem,
    const bool incrementColumnDir
)
{
    // Sanity checks
    {
        const direction blockMatrixSize = pTraits<blockType>::nComponents;
        const direction blockMatrixASize = pTraits<Type>::nComponents;

        if (blockMatrixSize > blockMatrixASize)
        {
            FatalErrorIn
            (
                "void fvBlockMatrix<Type>::insertBlock\n"
                "(\n"
                "    const direction dirI,\n"
                "    const direction dirJ,\n"
                "    BlockLduSystem<blockType, fieldType>& blockSystem,\n"
                "    const bool incrementColumnDir\n"
                ")"
            )   << "Trying to insert a block matrix from BlockLduSystem into "
                << "smaller one from fvBlockMatrix."
                << abort(FatalError);
        }

        if (dirI == dirJ)
        {
            FatalErrorIn
            (
                "void fvBlockMatrix<Type>::insertBlock\n"
                "(\n"
                "    const direction dirI,\n"
                "    const direction dirJ,\n"
                "    BlockLduSystem<blockType, fieldType>& blockSystem,\n"
                "    const bool incrementColumnDir\n"
                ")"
            )   << "Trying to insert coupling in the position where equation "
                << "should be, since dirI = dirJ. Try using insertEquation "
                << "member function."
                << abort(FatalError);
        }
    }

    const direction nCmpts = pTraits<blockType>::nComponents;
    direction localDirI = dirI;
    direction localDirJ = dirJ;

    // Get references to ldu fields of blockMatrix always as linear
    const typename CoeffField<blockType>::linearTypeField& bmd =
        blockSystem.diag().asLinear();
    const typename CoeffField<blockType>::linearTypeField& bmu =
        blockSystem.upper().asLinear();
    const typename CoeffField<blockType>::linearTypeField& bml =
        blockSystem.lower().asLinear();

    // Get references to ldu fields of this block matrix always as square
    typename CoeffField<Type>::squareTypeField& blockDiag =
         this->diag().asSquare();
    typename CoeffField<Type>::squareTypeField& blockUpper =
         this->upper().asSquare();
    typename CoeffField<Type>::squareTypeField& blockLower =
         this->lower().asSquare();

    // Insert blockMatrix that represents coupling into larger system matrix
    for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
    {
        forAll(bmd, cellI)
        {
            blockDiag[cellI](localDirI, localDirJ) +=
                bmd[cellI].component(cmptI);
        }

        forAll(bmu, faceI)
        {
            blockUpper[faceI](localDirI, localDirJ) +=
                bmu[faceI].component(cmptI);
            blockLower[faceI](localDirI, localDirJ) +=
                bml[faceI].component(cmptI);
        }

        if (incrementColumnDir)
        {
            localDirI++;
        }
        else
        {
            localDirJ++;
        }
    }
}


template<class Type>
template<class blockType, class fieldType>
void Foam::fvBlockMatrix<Type>::insertBoundaryContributions
(
    const direction dirI,
    const direction dirJ,
    const BlockLduSystem<blockType, fieldType>& blockSystem,
    const bool incrementColumnDir
)
{
    // Need to get reference to fvMesh instead of lduMesh
    const fvBoundaryMesh& bmesh =
        refCast<const fvMesh>(this->mesh()).boundary();

    const direction nCmpts = pTraits<blockType>::nComponents;
    const direction nSrcCmpts = pTraits<fieldType>::nComponents;
    direction localDirI = dirI;
    direction localDirJ = dirJ;

    const Field<fieldType>& source = blockSystem.source();

    // Get reference to this source field of block system
    Field<Type>& b = this->source();

    // Insert source from block system to this system's rhs
    for (direction cmptI = 0; cmptI < nSrcCmpts; cmptI++)
    {
        scalarField sourceCmpt(source.component(cmptI));

        forAll(b, cellI)
        {
            b[cellI](localDirI) += sourceCmpt[cellI];
        }

        if (incrementColumnDir)
        {
            localDirI++;
        }
        else
        {
            localDirJ++;
        }
    }

    // Reset local directions for coupling contributions
    localDirI = dirI;
    localDirJ = dirJ;

    // Insert coupling contributions into block matrix
    forAll(bmesh, patchI)
    {
        if (bmesh[patchI].coupled())
        {
            typename CoeffField<Type>::squareTypeField& pcoupleUpper =
                this->coupleUpper()[patchI].asSquare();
            typename CoeffField<Type>::squareTypeField& pcoupleLower =
                this->coupleLower()[patchI].asSquare();

            const typename CoeffField<blockType>::linearTypeField& bmcu =
                blockSystem.coupleUpper()[patchI].asLinear();
            const typename CoeffField<blockType>::linearTypeField& bmcl =
                blockSystem.coupleLower()[patchI].asLinear();

            for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
            {
                forAll(bmcu, faceI)
                {
                    pcoupleUpper[faceI](localDirI, localDirJ) +=
                          bmcu[faceI].component(cmptI);
                    pcoupleLower[faceI](localDirI, localDirJ) +=
                          bmcl[faceI].component(cmptI);
                }

                if (incrementColumnDir)
                {
                    localDirI++;
                }
                else
                {
                    localDirJ++;
                }
            }

            // Reset local directions for other patches
            localDirI = dirI;
            localDirJ = dirJ;
        }
    }
}


template<class Type>
void fvBlockMatrix<Type>::insertCouplingDiag
(
    const direction dirI,
    const direction dirJ,
    const scalarField& coeffIJ
)
{
    // Get reference to block diagonal of the block system
    typename CoeffField<Type>::squareTypeField& blockDiag =
        this->diag().asSquare();

    // Set off-diagonal coefficient
    forAll(coeffIJ, cellI)
    {
        blockDiag[cellI](dirI, dirJ) += coeffIJ[cellI];
    }

    // Source compensation is done in function updateSourceCoupling()
    // after all coupling terms are added.  HJ, 27/Apr/2015
}


template<class Type>
void fvBlockMatrix<Type>::insertCouplingUpperLower
(
    const direction dirI,
    const direction dirJ,
    const fvScalarMatrix& matrix
)
{
    if (matrix.diagonal())
    {
        // Matrix for insertion is diagonal-only: nothing to do
        return;
    }

    if (matrix.symmetric() && this->symmetric())
    {
        Info<< "Both fvScalarMatrix and block matrix are symmetric: " << nl
            << "inserting only upper triangle"
            << endl;
    }
    else
    {
        // Either scalar or block matrix is asymmetric: insert lower triangle
        const scalarField& lower = matrix.lower();

        typename CoeffField<Type>::squareTypeField& blockLower =
            this->lower().asSquare();

        forAll (lower, cellI)
        {
            blockLower[cellI](dirI, dirJ) = lower[cellI];
        }
    }

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();

        typename CoeffField<Type>::squareTypeField& blockUpper =
            this->upper().asSquare();

        forAll(upper, cellI)
        {
            blockUpper[cellI](dirI, dirJ) = upper[cellI];
        }
    }
    else
    {
        FatalErrorIn
        (
            "void fvBlockMatrix<Type>::insertCouplingUpperLower\n"
            "(\n"
            "    const direction dirI,\n"
            "    const direction dirJ,\n"
            "    const fvScalarMatrix& m\n"
            ")"
        )   << "Error in matrix insertion: problem with block structure."
            << abort(FatalError);
    }

    // Insert block interface fields
    forAll(this->interfaces(), patchI)
    {
        if (this->interfaces().set(patchI))
        {
            // Couple upper and lower
            const scalarField& cUpper = matrix.boundaryCoeffs()[patchI];
            const scalarField& cLower = matrix.internalCoeffs()[patchI];

            typename CoeffField<Type>::squareTypeField& blockUpper =
                this->coupleUpper()[patchI].asSquare();

            typename CoeffField<Type>::squareTypeField& blockLower =
                this->coupleLower()[patchI].asSquare();

            forAll(cUpper, faceI)
            {
                blockUpper[faceI](dirI, dirJ) = cUpper[faceI];
                blockLower[faceI](dirI, dirJ) = cLower[faceI];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
fvBlockMatrix<Type>::fvBlockMatrix
(
    GeometricField<Type, fvPatchField, volMesh>& psi
)
:
    BlockLduSystem<Type, Type>(psi.mesh()),
    psi_(psi)
{
    this->interfaces() = psi.boundaryField().blockInterfaces();
}


template<class Type>
fvBlockMatrix<Type>::fvBlockMatrix
(
    const fvBlockMatrix<Type>& bxs
)
:
    BlockLduSystem<Type, Type>(bxs),
    psi_(bxs.psi())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
template<class fieldType>
void fvBlockMatrix<Type>::retrieveSolution
(
    const direction dir,
    Field<fieldType>& xSingle
) const
{
    const direction nCmpts = pTraits<fieldType>::nComponents;
    direction localDir = dir;

    const Field<Type> psiIn = psi_.internalField();

    for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
    {
        scalarField xSingleCurr(xSingle.component(cmptI));

        forAll (xSingleCurr, cellI)
        {
            xSingleCurr[cellI] = psiIn[cellI](localDir);
        }

        xSingle.replace(cmptI, xSingleCurr);

        localDir++;
    }
}


template<class Type>
template<class matrixType>
void fvBlockMatrix<Type>::insertEquation
(
    const direction dir,
    fvMatrix<matrixType>& matrix
)
{
    insertSolutionVector(dir, matrix.psi().internalField());
    insertDiagSource(dir, matrix);
    insertUpperLower(dir, matrix);
    updateCouplingCoeffs(dir, matrix);
}


template<class Type>
template<class blockType, class fieldType>
void fvBlockMatrix<Type>::insertBlockCoupling
(
    const direction dirI,
    const direction dirJ,
    const BlockLduSystem<blockType, fieldType>& blockSystem,
    const bool incrementColumnDir
)
{
    insertBlock(dirI, dirJ, blockSystem, incrementColumnDir);
    insertBoundaryContributions(dirI, dirJ, blockSystem, incrementColumnDir);
}


template<class Type>
void Foam::fvBlockMatrix<Type>::insertEquationCoupling
(
    const direction dirI,
    const direction dirJ,
    const scalarField& coeffIJ
)
{
    // Multiply coefficients by volume
    scalarField coeffIJVol = coeffIJ*psi_.mesh().V();

    insertCouplingDiag(dirI, dirJ, coeffIJVol);
}


template<class Type>
void fvBlockMatrix<Type>::blockAdd
(
    const direction dir,
    const scalarField& xSingle,
    Field<Type>& blockX
)
{
    forAll(xSingle, cellI)
    {
        blockX[cellI](dir) += xSingle[cellI];
    }
}


template<class Type>
void Foam::fvBlockMatrix<Type>::updateSourceCoupling()
{
    // Eliminate off-diagonal block coefficients from the square diagonal
    // With this change, the segregated matrix can be assembled with complete
    // source terms and linearisation can be provided independently.
    // Once the linearisation coefficients are set (off-diagonal entries
    // in the square block matrix, they are multiplied by the current value
    // of the field and subtracted from the source term

    if (this->diag().activeType() == blockCoeffBase::SQUARE)
    {
        typename CoeffField<Type>::squareTypeField& blockDiag =
            this->diag().asSquare();

        typename CoeffField<Type>::linearTypeField lf(blockDiag.size());
        typename CoeffField<Type>::squareTypeField sf(blockDiag.size());

        // Expand and contract

        // Take out the diagonal entries from the square coefficient
        contractLinear(lf, blockDiag);

        // Expand the diagonal for full square, with zeroes in the off-diagonal
        expandLinear(sf, lf);

        // Subtract from the source the difference between the full block
        // diagonal and the diagonal terms only
        // Sign is the same as in the derivative
        this->source() += (blockDiag - sf) & psi_.internalField();
    }
}


template<class Type>
BlockSolverPerformance<Type> fvBlockMatrix<Type>::solve
(
    const dictionary& solverControls
)
{
    // Solver call
    BlockSolverPerformance<Type> solverPerf =
        BlockLduSolver<Type>::New
        (
            psi_.name(),
            *this,
            solverControls
        )->solve(psi_.internalField(), this->source());

    // Print performance
    solverPerf.print();

    return solverPerf;
}


template<class Type>
BlockSolverPerformance<Type> fvBlockMatrix<Type>::solve()
{
    return solve(psi_.mesh().solutionDict().solverDict(psi_.name()));
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Ostream& operator<<
(
    Ostream& os,
    const fvBlockMatrix<Type>& bxs
)
{
    os  << static_cast<const BlockLduSystem<Type, Type>&>(bxs) << nl
        << bxs.psi() << endl;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
