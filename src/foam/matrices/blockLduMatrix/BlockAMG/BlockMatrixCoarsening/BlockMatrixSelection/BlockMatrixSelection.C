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
    BlockMatrixSelection

Description
    Selective AMG policy for block matrices

Author
    Tessa Uroic, FMENA

\*----------------------------------------------------------------------------*/

#include "BlockMatrixSelection.H"
#include "coeffFields.H"
#include "BlockAMGInterfaceField.H"
#include "coarseBlockAMGLevel.H"
#include "PriorityList.H"
#include "crMatrix.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Factor for defining strong negative-coupling of variables
template<class Type>
const Foam::debug::tolerancesSwitch
Foam::BlockMatrixSelection<Type>::epsilon_
(
    "samgWeightFactor",
    0.25
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::BlockMatrixSelection<Type>::calcCoarsening()
{
//------------------------------------------------------------------------------
//             MATRIX DATA: ADDRESSING, COEFFICIENTS, COEFF NORMS
//------------------------------------------------------------------------------
    Info<< "Start equation selection" << endl;
    // Get addressing
    const unallocLabelList& rowStart = matrix_.lduAddr().ownerStartAddr();
    const unallocLabelList& losortAddr = matrix_.lduAddr().losortAddr();
    const unallocLabelList& losortStart = matrix_.lduAddr().losortStartAddr();
    const unallocLabelList& column = matrix_.lduAddr().upperAddr();
    const unallocLabelList& row = matrix_.lduAddr().lowerAddr();
    const label matrixSize = matrix_.lduAddr().size();

    // Note: not taking norm magnitudes.  HJ, 28/Feb/2017

    // Calculate norm for diagonal coefficients
    scalarField normDiag(matrixSize);
    normPtr_->normalize(normDiag, matrix_.diag());

    // Note: this needs to be untangled for symmetric matrices.
    // If the matrix is symmetric and you ask for lower, it will be manufactured
    // for you and will double the memory.  Therefore, all loops dealing with
    // lower need to be under (if matrix.assymetric()) {...} protection.
    // Please get the code to work first and then refactor
    // HJ, 18/Feb/2017

    // Calculate norm for upper triangle coeffs (magUpper)
    scalarField normUpper(row.size());
    normPtr_->normalize(normUpper, matrix_.upper());

    // Calculate norm for lower triangle coeffs (magLower)
    scalarField normLower(column.size());
    normPtr_->normalize(normLower, matrix_.lower());

    // Calculate norm magnitudes
    scalarField magNormUpper = mag(normUpper);
    scalarField magNormLower = mag(normLower);

//------------------------------------------------------------------------------
//        CRITERIA FOR COARSENING: FIND STRONG CONNECTIONS FOR EACH ROW
//------------------------------------------------------------------------------

    // Find the largest norm in the row - use it as a criterion for strong
    // connections
    // Compare the elements in each row to the largest coefficient stored
    // in largestNorm[i] multiplied by the constant epsilon and declare
    // all larger elements to be influences of i
    // NOTE: This approach may only be valid for pressure as the direction
    // because it has coeffs on the off-diagonal with sign opposite to diagonal

    // Create largest norm.  It will be multiplied by epsilon later
    scalarField epsLargestNorm(matrixSize, 0);

    // Count strong elements in each row - for rowStart addressing
    // Note: checking magnitude of coeff, which accounts for both strong
    // positive and negative coefficients.  HJ, 28/Feb/2017
    labelList strongCoeffCounter(matrixSize, 0);

    register label j;

    for (register label i = 0; i < matrixSize; i++)
    {
        for (register label k = rowStart[i]; k < rowStart[i+1]; k++)
        {
            // Upper triangle
            if (magNormUpper[k] > epsLargestNorm[i])
            {
                epsLargestNorm[i] = magNormUpper[k];
            }

            // Lower triangle
            j = column[k];

            if (magNormLower[k] > epsLargestNorm[j])
            {
                epsLargestNorm[j] = magNormLower[k];
            }
        }
    }

    // Multiply largest norm by epsilon.  This is now it is used below
    epsLargestNorm *= epsilon_();

    // Count strong elements in the matrix to create addressing
    label count = 0;

    for (register label i = 0; i < matrixSize; i++)
    {
        for (register label k = rowStart[i]; k < rowStart [i+1]; k++)
        {
            // Column of coefficient k in row i
            j = column[k];

            // Do the upper triangle (rowwise)
            if (magNormUpper[k] > epsLargestNorm[i])
            {
                    strongCoeffCounter[i]++;
                    count++;
            }

            // Do the lower triangle (columnwise)
            if (magNormLower[k] > epsLargestNorm[j])
            {
                strongCoeffCounter[j]++;
                count++;
            }
        }
    }


    // Create a crMatrix that will store all of the strong elements of each row
    // (some will become FINE and some COARSE - and this will determine the
    // direct and standard interpolation procedures)

    crMatrix strong(matrixSize, matrixSize, strongCoeffCounter);

    // Set addressing for the matrix:
    // stongCol and strongElement are arrays needed to create a compressed row
    // matrix
    const labelList& strongRowStart = strong.crAddr().row();
    labelList& strongCol = strong.col();
    scalarField& strongCoeff = strong.coeffs();

    // Counter for counting the strong elements in a row
    labelList counter(matrixSize, 0);

    // Counter for strong elements in a column
    labelList influence(matrixSize, 0);

    for (register label i = 0; i < matrixSize; i++)
    {
        for (register label k = rowStart[i]; k < rowStart [i+1]; k++)
        {
            j = column[k];

            // Check elements in upper triangle
            if (magNormUpper[k] > epsLargestNorm[i])
            {
                // Store the strong elements into crMatrix, use counter
                // to count the number of negative strong elements for
                // each row i
                strongCol[strongRowStart[i] + counter[i]] = j;
                strongCoeff[strongRowStart[i] + counter[i]] = normUpper[k];
                counter[i]++;
                influence[j]++;
            }

            // Check elements in lower triangle
            if (magNormLower[k] > epsLargestNorm[j])
            {
                strongCol[strongRowStart[j] + counter[j]] = i;
                strongCoeff[strongRowStart[j] + counter[j]] = normLower[k];
                counter[j]++;
                influence[i]++;
            }
        }
    }

//------------------------------------------------------------------------------
//           COARSENING: SORT EQUATIONS INTO COARSE AND FINE SUBSETS
//------------------------------------------------------------------------------

    // Transpose the compressed row matrix to use for coarsening
    crAddressing Taddr = strong.crAddr().T();
    const labelList& tRow = Taddr.row();
    const labelList& tCol = Taddr.col();

    // Label the equations COARSE and FINE based on the number of
    // influences.
    // In order to do that a priority list can be used to determine the
    // equation  with the (currently) largest number of influences
    // (weight of the element)

    // Mark equations
    rowLabel_.setSize(matrixSize, UNDECIDED);

    PriorityList<label> equationWeight(matrixSize);
    for (label i = 0; i < matrixSize; i++)
    {
        // Set weights for each equation (weight == number of strong
        // connections in column!)
        equationWeight.set(i, influence[i]);

        // Label rows without connections as FINE
        if (counter[i] == 0)
        {
            rowLabel_[i] = FINE;
        }
    }

    while (!equationWeight.empty())
    {
        // removeHead = return index with largest weight and remove
        label topElement = equationWeight.removeHead();

        if (rowLabel_[topElement] == UNDECIDED)
        {
            rowLabel_[topElement] = COARSE;
            nCoarseEqns_++;

            // UPPER TRIANGLE - row-wise
            for
            (
                label i = rowStart[topElement];
                i < rowStart[topElement + 1];
                i++
            )
            {
                label neighbour = column[i];

                if (rowLabel_[neighbour] == UNDECIDED)
                {
                    equationWeight.updateWeight(neighbour, -1);
                }
            }
            // LOWER TRIANGLE - row-wise (losort)
            for
            (
                label i = losortStart[topElement];
                i < losortStart[topElement + 1];
                i++
            )
            {
                label index = losortAddr[i];
                label neighbour = row[index];

                if (rowLabel_[neighbour] == UNDECIDED)
                {
                    equationWeight.updateWeight(neighbour, -1);
                }
            }

            for
            (
                label i = tRow[topElement];
                i < tRow[topElement + 1];
                i++
            )
            {
                label neighbour = tCol[i];

                if (rowLabel_[neighbour] == UNDECIDED)
                {
                    rowLabel_[neighbour] = FINE;

                    for
                    (
                        label j  = strongRowStart[neighbour];
                        j < strongRowStart[neighbour + 1];
                        j++
                    )
                    {
                        label acquaintance = strongCol[j];
                        if (rowLabel_[acquaintance] == UNDECIDED)
                        {
                            equationWeight.updateWeight(acquaintance, 2);
                        }
                    }
                }
            }
        }
    }

//------------------------------------------------------------------------------
//                     CALCULATE COARSE MATRIX ADDRESSING
//------------------------------------------------------------------------------

    // Calculate the number of fine rows with index smaller than i and
    // save it - it will be used for prolongation matrix addresing
    labelField fineRowsCount(matrixSize, 0);

    for (label i = 0; i < matrixSize; i++)
    {
        if (rowLabel_[i] == FINE)
        {
            if (i == 0)
            {
                fineRowsCount[i] = 1;
            }
            else
            {
                fineRowsCount[i] = fineRowsCount[i - 1] + 1;
            }
        }
        else
        {
            if (i == 0)
            {
                continue;
            }
            else
            {
                fineRowsCount[i] = fineRowsCount[i - 1];
            }
        }

    }

//------------------------------------------------------------------------------
//              CALCULATING CONTRIBUTIONS TO THE SCALING FACTOR
//------------------------------------------------------------------------------

    // Sum of positive elements in row (sign equal to diagonal)
    scalarField positiveElemSum(matrixSize, 0);

    // Sum of negative elements in row
    scalarField negativeElemSum(matrixSize, 0);

    // Sum of negative COARSE elements in row
    scalarField negCoarseElemSum(matrixSize, 0);

    // Sum of positive COARSE elements in row
    scalarField posCoarseElemSum(matrixSize, 0);

    // Adding positive and negative coeffs in each row (numerator of
    // scaling factor)
    for (label i = 0; i < matrixSize; i++)
    {
        // Only need neighbours for interpolation of FINE rows
        if (rowLabel_[i] == FINE)
        {
            // Extract diagonal component to compare the sign with the
            // off-diagonal coeffs - search for negative and positive
            // connections
            scalar diagComponent = normDiag[i];

            // Upper triangle - use rowStart
            for (label k = rowStart[i]; k < rowStart[i + 1]; k++)
            {
                // POSITIVE contribution
                if (sign(normUpper[k]) == sign(diagComponent))
                {
                    positiveElemSum[i] += normUpper[k];
                }
                // NEGATIVE contribution
                else
                {
                    negativeElemSum[i] += normUpper[k];
                }
            }

            // Lower triangle - use losort to go row-wise (Why? This
            // is a pain)
            // Because I have to compare the sign to the
            // corresponding diagonal element in row i.
            for (label k = losortStart[i]; k < losortStart[i + 1]; k++)
            {
                label index = losortAddr[k];

                // POSITIVE contribution
                if (sign(normLower[index]) == sign(diagComponent))
                {
                    positiveElemSum[i] += normLower[index];
                }
                // NEGATIVE contribution
                else
                {
                    negativeElemSum[i] += normLower[index];
                }
            }
        }
    }

    // Adding only strong coarse contributions for each row (denominator
    // of scaling factor)

    for (label i = 0; i < matrixSize; i++)
    {
        if (rowLabel_[i] == FINE)
        {
            scalar diagSign = sign(normDiag[i]);
            for
            (
                label index = strongRowStart[i];
                index < strongRowStart[i + 1];
                index++
            )
            {
                // Get column of strong coeff
                label j = strongCol[index];

                // Sum positive connections
                if
                (
                    sign(strongCoeff[index]) == diagSign &&
                    rowLabel_[j] == COARSE
                )
                {
                    posCoarseElemSum[i] += strongCoeff[index];
                }
                else if
                (
                    rowLabel_[j] == COARSE
                )
                {
                    negCoarseElemSum[i] += strongCoeff[index];
                }
            }
        }
    }

//------------------------------------------------------------------------------
//                    CALCULATING PROLONGATION WEIGHTS
//------------------------------------------------------------------------------

    // Fields for assembling prolongation matrix
    // We don't know in advance the number of coeffs in the
    // prolongation, but per row, there cannot be more than the number
    // of coarse equations
    labelList prolongationRow(matrixSize + 1);
    DynamicList<label> prolongationCol;
    DynamicList<scalar> prolongationCoeff;

    // Initialize starting row
    prolongationRow[0] = 0;

    for (label i = 0; i < matrixSize; i++)
    {
        label rowCount = prolongationRow[i];

        if (rowLabel_[i] == FINE)
        {
            scalar diagonalCoeff = normDiag[i];

            // Check whether this Eqn has COARSE negative strong
            // contributions!

            if (negCoarseElemSum[i] == 0)
            {
#   ifdef FULLDEBUG

                // Positive connections are used for interpolation! Diagonal
                // dominance is not conserved!
                Info << "Interpolation from positive neighbours!" << nl
                     << "Equation " << i << endl;

#   endif

                // Row has no negative COARSE contributions.
                // Interpolating from positive contributions!
                for
                (
                    label k = strongRowStart[i];
                    k < strongRowStart[i + 1];
                    k++
                )
                {
                    label j = strongCol[k];

                    if (rowLabel_[j] == COARSE)
                    {
                        prolongationCoeff.append
                        (
                            (-positiveElemSum[i])
                           /posCoarseElemSum[i]
                           *strongCoeff[k]
                           /(diagonalCoeff)
                        );

                        prolongationCol.append(j - fineRowsCount[j]);
                        rowCount++;

                    }
                }
            }
            else
            {
                // Row has negative COARSE contributions. Use only those
                // for interpolation!
                for
                (
                    label k = strongRowStart[i];
                    k < strongRowStart [i + 1];
                    k++
                )
                {
                    label j = strongCol[k];

                    if
                    (
                        rowLabel_[j] == COARSE &&
                        sign(strongCoeff[k]) != sign(diagonalCoeff)
                    )
                    {
                        prolongationCoeff.append
                        (
                            (-negativeElemSum[i])
                           /negCoarseElemSum[i]
                           *strongCoeff[k]
                           /(diagonalCoeff + positiveElemSum[i])
                        );

                        prolongationCol.append(j - fineRowsCount[j]);
                        rowCount++;
                    }
                }

            }
        }
        else // The row is COARSE - injection
        {
            prolongationCoeff.append(1);
            prolongationCol.append(i - fineRowsCount[i]);
            rowCount++;
        }

        // Go to next row in the next loop
        prolongationRow[i + 1] = rowCount;
    }

    // Resize column array and coefficients array
    prolongationCoeff.setSize(prolongationRow[matrixSize]);
    prolongationCol.setSize(prolongationRow[matrixSize]);

    // Assemble prolongation matrix
    Pptr_ = new crMatrix
    (
        matrixSize,
        nCoarseEqns_,
        prolongationRow,
        prolongationCol
    );
    crMatrix& prolongation = *Pptr_;

    prolongation.coeffs().transfer(prolongationCoeff);

    // Ranking for smoothing sweeps - first all coarse equations,
    // then all fine
    labelField eqnRank_(matrixSize, 0);
    label rank = 0;

    // COARSE
    for (label i = 0; i < matrixSize; i++)
    {
        if (rowLabel_[i] == COARSE)
        {
            eqnRank_[rank] = i;
            rank++;
        }
    }

    // FINE
    for (label i = 0; i < matrixSize; i++)
    {
        if (rowLabel_[i] == FINE)
        {
            eqnRank_[rank] = i;
            rank++;
        }
    }

    // If the number of coarse equations is less than minimum and
    // if the matrix has reduced in size by at least 1/3, coarsen
    if
    (
        nCoarseEqns_ > this->minCoarseEqns()
     && 3*nCoarseEqns_ <= 2*matrixSize
    )
    {
        coarsen_ = true;
    }

    reduce(coarsen_, andOp<bool>());

    if (coarsen_)
    {
        // Coarsening OK, make restriction matrix

        // Create restriction by transposing prolongation
        Rptr_ = new crMatrix(Pptr_->T());
    }
    else
    {
        // Coarsening did not succeed.  Delete Pptr
        deleteDemandDrivenData(Pptr_);
    }

    Info<< "End equation selection" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockMatrixSelection<Type>::BlockMatrixSelection
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict,
    const label groupSize,
    const label minCoarseEqns
)
:
    BlockMatrixCoarsening<Type>(matrix, dict, groupSize, minCoarseEqns),
    matrix_(matrix),
    normPtr_(BlockCoeffNorm<Type>::New(dict)),
    nCoarseEqns_(0),
    coarsen_(false),
    Pptr_(NULL),
    Rptr_(NULL),
    rowLabel_()
{
    calcCoarsening();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockMatrixSelection<Type>::~BlockMatrixSelection()
{
    deleteDemandDrivenData(Rptr_);
    deleteDemandDrivenData(Pptr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockAMGLevel<Type> >
Foam::BlockMatrixSelection<Type>::restrictMatrix() const
{
    Info<< "Start matrix restriction" << endl;
    if (!coarsen_)
    {
        FatalErrorIn("autoPtr<amgMatrix> samgPolicy::restrictMatrix() const")
            << "Requesting coarse matrix when it cannot be created"
            << abort(FatalError);
    }

    // Get references to restriction and prolongation matrix

    const crMatrix& R = *Rptr_;
    const crMatrix& P = *Pptr_;

    // Get connectivity
    const crAddressing& crR = R.crAddr();
    const crAddressing& crP = P.crAddr();

#   ifdef FULLDEBUG
    // Check sized chain rule
    const label nEqns = matrix_.lduAddr().size();

    if
    (
        crR.nCols() != nEqns
     || crP.nRows() != nEqns
     || crR.nRows() != nCoarseEqns_
    )
    {
        FatalErrorIn
        (
            "autoPtr<Foam::BlockAMGLevel<Type> >"
            "BlockMatrixSelection<Type>::restrictMatrix() const"
        )   << "Incompatible matrices for triple product: "
            << "R( " << crR.nRows() << " ," << crR.nCols() << ") "
            << "A( " << nEqns << " ," << nEqns << ") "
            << "P( " << crP.nRows() << " ," << crP.nCols() << ") "
            << abort(FatalError);
    }
#   endif


//------------------------------------------------------------------------------
//                                GET ADDRESSING
//------------------------------------------------------------------------------

    // Restriction addressing
    const labelList& rowStartR = crR.row();
    const labelList& colR = crR.col();

    // Matrix A addressing
    const unallocLabelList& rowStartA = matrix_.lduAddr().ownerStartAddr();
    const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();
    const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
    const unallocLabelList& losortAddr = matrix_.lduAddr().losortAddr();
    const unallocLabelList& losortStart = matrix_.lduAddr().losortStartAddr();

    // Prolongation addressing
    const labelList& rowStartP = crP.row();
    const labelList& colP = crP.col();

    // coeffLabel is used for checking if a contribution in that address already
    // exists
    labelList coeffLabel(nCoarseEqns_, -1);
    label nCoarseCoeffs = 0;

    // Array for creating rowStart addressing for upper triangle
    labelList nUpperCoeffs(nCoarseEqns_, 0);

    // Create coarse addressing - owner and neighbour pairs
    // HJ: give a better guess of size of coarse off-diagonal
    DynamicList<label> coarseOwner(nCoarseEqns_);
    DynamicList<label> coarseNeighbour(nCoarseEqns_);

//------------------------------------------------------------------------------
//     COUNT COARSE COEFFICIENTS IN TRIPLE PRODUCT AND CREATE ADDRESSING
//------------------------------------------------------------------------------

    // This loop will be used to count the number of coeffs appearing in the
    // COARSE matrix and creating the owner-neighbour addressing of the coarse
    // matrix.

    // Dynamic Lists are used to create coarseOwner and coarseNeighbour arrays
    // in which owner and neighbour are stored for the upper tringle. Note: the
    // neighbour array should be sorted before constructing the coarseMatrix.

    // TRIPLE PRODUCT:
    //     R*A multiplication
    //         Multiply matrix A coeff with R coeff if rowA == columnR
    //         The address of the resulting coeff is (rowR, columnA)
    //     RA*P multiplication
    //         Multiply P coeff with RA coeff if rowP == columnRA
    //         The address of the resulting coeff is (rowRA, columnP), that is
    //         the address in the COARSE matrix is (rowR, columnP)

    // NOTE: Letters i and j are used to denote the row and the column of the
    // matrix, respectively. In addition, row of matrix A is denoted ia, column
    // ja, row of prolongation is ip, etc.

    // Loop through rows of R
    for (label ir = 0; ir < nCoarseEqns_; ir++)
    {
        // Compressed row format - get indices of coeffsR in row ir
        for (label indexR = rowStartR[ir]; indexR < rowStartR[ir + 1]; indexR++)
        {
            // Column of R coeff
            const label jr = colR[indexR];

            // LOWER TRIANGLE
            for
            (
                label indexA = losortStart[jr];
                indexA < losortStart[jr + 1];
                indexA++
            )
            {
                // Get face index of coeff in A
                const label faceA = losortAddr[indexA];

                // Column of coeff in A
                const label ja = lowerAddr[faceA];

                // Go into the corresponding row of prolongation to find
                // contributions
                for
                (
                    label indexP = rowStartP[ja];
                    indexP < rowStartP[ja + 1];
                    indexP++
                )
                {
                    // Column of coeff in P
                    const label jp = colP[indexP];

                    // Check the address of the coefficient in COARSE matrix
                    // (ir = row, jp = column)
                    if (ir > jp)
                    {
                        // This coeff belongs to the lower triangle, and
                        // corresponding to lduAddressing, its address won't be
                        // stored
                        continue;
                    }
                    else if (ir != jp)
                    {
                        // Does any contribution to this address already exist?
                        // Label the column with the row index in which we found
                        // the coeff. If there is another one with this address,
                        // the column is already labeled with the row index.
                        if (coeffLabel[jp] != ir)
                        {
                            coeffLabel[jp] = ir;

                            coarseOwner.append(ir);
                            coarseNeighbour.append(jp);

                            nCoarseCoeffs++;

                            // Count upper triangle coeff for this row
                            nUpperCoeffs[ir]++;
                        }
                    }
                }
            }

            // UPPER TRIANGLE
            for
            (
                label faceA  = rowStartA[jr];
                faceA < rowStartA[jr + 1];
                faceA++
            )
            {
                // Get column of coeff in A
                const label ja = upperAddr[faceA];

                // Go into the corresponding row of prolongation to find
                // contributions
                for
                (
                    label indexP = rowStartP[ja];
                    indexP < rowStartP[ja + 1];
                    indexP++
                )
                {
                    // Get column of coefficient in P
                    const label jp = colP[indexP];

                    // Check the address of the coefficient in COARSE matrix
                    // (ir = row, jp = column)
                    if (ir > jp)
                    {
                        // This coeff belongs to the lower triangle, and
                        // corresponding to lduAddressing, its address won't be
                        // stored
                        continue;
                    }
                    // Count the coeff if it is an off-diagonal coeff in the
                    // upper triangle
                    else if (ir != jp)
                    {
                        // Does any contribution to this address already exist?
                        // We will label the column with the row index in which
                        // we found the coeff. If there is another one, the
                        // column is already labeled with the row index.
                        if (coeffLabel[jp] != ir)
                        {
                            coeffLabel[jp] = ir;

                            coarseOwner.append(ir);
                            coarseNeighbour.append(jp);

                            nCoarseCoeffs++;

                            // Count upper triangle coeff for this row
                            nUpperCoeffs[ir]++;
                        }
                    }
                }
            }

            // DIAGONAL
            for
            (
                label indexP = rowStartP[ir];
                indexP < rowStartP[ir + 1];
                indexP++
            )
            {
                // Column of coefficient in P
                const label jp = colP[indexP];

                // Check the address of the coefficient
                // (ir = row, jp = column)
                if (ir > jp)
                {
                    // This coeff belongs to the lower triangle, and
                    // corresponding to lduAddressing, its address won't be
                    // stored
                    continue;
                }
                // Count the coeff if it is an off-diagonal coeff in the
                // upper triangle
                else if (ir != jp)
                {
                    // Does any contribution to this address already exist?
                    // We will label the column with the row index in which
                    // we found the coeff. If there is another one, the
                    // column is already labeled with the row index.
                    if (coeffLabel[jp] != ir)
                    {
                        coeffLabel[jp] = ir;

                        coarseOwner.append(ir);
                        coarseNeighbour.append(jp);

                        nCoarseCoeffs++;

                        // Count upper triangle coeff for this row
                        nUpperCoeffs[ir]++;
                    }
                }
            }
        }
    }

//------------------------------------------------------------------------------
//                      CREATE COARSE MATRIX ADDRESSING
//------------------------------------------------------------------------------

    // Create rowStartAddressing for upperCoeffs
    labelList coarseRowStart(nCoarseEqns_ + 1, 0);
    for (label i = 1; i <= nCoarseEqns_; i++)
    {
        coarseRowStart[i] = coarseRowStart[i - 1] + nUpperCoeffs[i - 1];
    }

    // Set size of dynamic list
    coarseOwner.setSize(nCoarseCoeffs);
    coarseNeighbour.setSize(nCoarseCoeffs);

    // Sorting the coarse neighbour list to go in ascending order
    for (label i = 0; i < nCoarseEqns_; i++)
    {
        SortableList<label> sortNeighbour(nUpperCoeffs[i]);
        label count = 0;

        for (label k = coarseRowStart[i]; k < coarseRowStart[i + 1]; k++)
        {
            sortNeighbour[count] = coarseNeighbour[k];
            count++;
        }

        // Sort neighbour list
        sortNeighbour.sort();

        // Reset count
        count = 0;

        // Copy sorted address indices into coarseNeighbour list
        for (label k = coarseRowStart[i]; k < coarseRowStart[i + 1]; k++)
        {
            coarseNeighbour[k] = sortNeighbour[count];
            count++;
        }
    }

    // Set the coarse ldu addressing onto the list
    autoPtr<lduPrimitiveMesh> coarseAddrPtr
    (
        new lduPrimitiveMesh
        (
            nCoarseEqns_,
            coarseOwner,
            coarseNeighbour,
            true
        )
    );

    const typename BlockLduInterfaceFieldPtrsList<Type>::Type&
        interfaceFields =
        const_cast<BlockLduMatrix<Type>&>(matrix_).interfaces();

    // Set the coarse interfaces and coefficients
    lduInterfacePtrsList coarseInterfaces(interfaceFields.size());

    labelListList coarseInterfaceAddr(interfaceFields.size());

    // Initialise transfer of restrict addressing on the interface
    // HJ, reconsider blocking comms.  HJ, 9/Jun/2016
    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            interfaceFields[intI].coupledInterface().initInternalFieldTransfer
            (
                Pstream::blocking,
                rowLabel_
            );
        }
    }

    // Store coefficients to avoid tangled communications
    // HJ, 1/Apr/2009
    FieldField<Field, label> fineInterfaceAddr(interfaceFields.size());

    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            const lduInterface& fineInterface =
                interfaceFields[intI].coupledInterface();

            fineInterfaceAddr.set
            (
                intI,
                new labelField
                (
                    fineInterface.internalFieldTransfer
                    (
                        Pstream::blocking,
                        rowLabel_
                    )
                )
            );
        }
    }

    // Create AMG interfaces
    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            const lduInterface& fineInterface =
                interfaceFields[intI].coupledInterface();

            coarseInterfaces.set
            (
                intI,
                AMGInterface::New
                (
                    coarseAddrPtr(),
                    coarseInterfaces,
                    fineInterface,
                    fineInterface.interfaceInternalField(rowLabel_),
                    fineInterfaceAddr[intI]
                ).ptr()
            );
        }
    }

    //HJ: Add interface fields
//------------------------------------------------------------------------------
//                            CREATE COARSE MATRIX
//------------------------------------------------------------------------------

    // Set the coarse level matrix
    autoPtr<BlockLduMatrix<Type> > coarseMatrixPtr
    (
        new BlockLduMatrix<Type>(coarseAddrPtr())
    );
    BlockLduMatrix<Type>& coarseMatrix = coarseMatrixPtr();

    typedef CoeffField<Type> TypeCoeffField;

    TypeCoeffField& coarseUpper = coarseMatrix.upper();
    TypeCoeffField& coarseDiag = coarseMatrix.diag();
    TypeCoeffField& coarseLower = coarseMatrix.lower();

    // Addresing of coarse matrix
    const unallocLabelList& coarseColumn = coarseMatrix.lduAddr().upperAddr();
    // const unallocLabelList& coarseRow = coarseMatrix.lduAddr().lowerAddr();


//------------------------------------------------------------------------------
//                           GET COEFFICIENTS
//------------------------------------------------------------------------------

    // Coefficients of matrix A
    const TypeCoeffField& diag = matrix_.diag();
    const TypeCoeffField& upper = matrix_.upper();
    const TypeCoeffField& lower = matrix_.lower();

    // Coefficients of restriction R
    const scalarField& coeffR = R.coeffs();

    // Coefficients of prolongation P
    const scalarField& coeffP = P.coeffs();

//------------------------------------------------------------------------------
//                TRIPLE PRODUCT = Restriction*A*Prolongation
//------------------------------------------------------------------------------

    // NOTE: Letters i and j are used to denote the row and the column of the
    // matrix, respectively. In addition, row of matrix A is denoted ia, column
    // ja, row of prolongation is ip, etc.

    if
    (
        diag.activeType() == blockCoeffBase::SQUARE
     && upper.activeType() == blockCoeffBase::SQUARE
    )
    {
        typedef typename TypeCoeffField::squareTypeField squareTypeField;
        typedef typename TypeCoeffField::squareType squareType;

        const squareTypeField& activeUpper = upper.asSquare();
        const squareTypeField& activeDiag = diag.asSquare();
        const squareTypeField& activeLower = lower.asSquare();

        squareTypeField& activeCoarseUpper = coarseUpper.asSquare();
        squareTypeField& activeCoarseDiag = coarseDiag.asSquare();
        squareTypeField& activeCoarseLower = coarseLower.asSquare();

        // Loop through rows of R
        for (label ir = 0; ir < nCoarseEqns_; ir++)
        {
            // Compressed row format, get indices of coeffsR in row ir
            for
            (
                label indexR = rowStartR[ir];
                indexR < rowStartR[ir + 1];
                indexR++
            )
            {
                // Column of coeff in R
                const label jr = colR[indexR];

                // FINE LOWER TRIANGLE
                for
                (
                    label indexA = losortStart[jr];
                    indexA < losortStart[jr + 1];
                    indexA++
                )
                {
                    // Get face index of coeff in A
                    const label faceA = losortAddr[indexA];

                    // Column index of coeff in A
                    const label ja = lowerAddr[faceA];

                    // Multiply coefficients of R and A
                    // Address is (rowR, columnA)
                    const squareType& ra = activeLower[faceA]*coeffR[indexR];

                    // Go into corresponding row of prolongation to find
                    // contributions
                    for
                    (
                        label indexP = rowStartP[ja];
                        indexP < rowStartP[ja + 1];
                        indexP++
                    )
                    {
                        // Column of coeff in P
                        const label jp = colP[indexP];

                        // Check the address of the coefficient
                        // (ir = row, jp = column)
                        if (ir > jp)
                        {
                            // Found lower COARSE triangle
                            // Find the corresponding owner-neighbour pair in
                            // the upper triangle
                            for
                            (
                                label m = coarseRowStart[jp];
                                m < coarseRowStart[jp + 1];
                                m++
                            )
                            {
                                if (ir == coarseColumn[m])
                                {
                                    activeCoarseLower[m] += ra*coeffP[indexP];
                                    break;
                                }
                            }
                        }
                        else if (ir == jp)
                        {
                            // Found COARSE diagonal coefficient
                            activeCoarseDiag[ir] += ra*coeffP[indexP];
                        }
                        else
                        {
                            // Found upper COARSE triangle
                            // Search for coefficient insertion.  HJ: this needs optimisation
                            for
                            (
                                label m = coarseRowStart[ir];
                                m < coarseRowStart[ir + 1];
                                m++
                            )
                            {
                                if (jp == coarseColumn[m])
                                {
                                    activeCoarseUpper[m] += ra*coeffP[indexP];
                                    break;
                                }
                            }
                        }
                    }
                }

                // FINE UPPER TRIANGLE
                for
                (
                    label faceA = rowStartA[jr];
                    faceA < rowStartA[jr + 1];
                    faceA++
                )
                {
                    const label ja = upperAddr[faceA];

                    const squareType& ra = coeffR[indexR]*activeUpper[faceA];

                    for
                    (
                        label indexP = rowStartP[ja];
                        indexP < rowStartP[ja + 1];
                        indexP++
                    )
                    {
                        const label jp = colP[indexP];

                        if (ir > jp)
                        {
                            // Found lower COARSE triangle
                            // Find owner-neighbour in the upper triangle and
                            // store
                            for
                            (
                                label m = coarseRowStart[jp];
                                m < coarseRowStart[jp + 1];
                                m++
                            )
                            {
                                if (ir == coarseColumn[m])
                                {
                                    activeCoarseLower[m] += ra*coeffP[indexP];
                                    break;
                                }
                            }
                        }
                        else if (ir == jp)
                        {
                            // Found COARSE diagonal
                            activeCoarseDiag[ir] += ra*coeffP[indexP];
                        }
                        else
                        {
                            // Found upper COARSE triangle
                            // Search for coefficient insertion.  HJ: this needs optimisation
                            for
                            (
                                label m = coarseRowStart[ir];
                                m < coarseRowStart[ir + 1];
                                m++
                            )
                            {
                                if (jp == coarseColumn[m])
                                {
                                    activeCoarseUpper[m] += ra*coeffP[indexP];
                                    break;
                                }
                            }
                        }
                    }
                }

                // FINE DIAGONAL
                const squareType& ra = coeffR[indexR]*activeDiag[jr];

                for
                (
                    label indexP= rowStartP[jr];
                    indexP < rowStartP[jr + 1];
                    indexP++
                )
                {
                    const label jp = colP[indexP];

                    if (ir > jp)
                    {
                        // Found lower COARSE triangle
                        // Search for coefficient insertion.  HJ: this needs optimisation
                        for
                        (
                            label m = coarseRowStart[jp];
                            m < coarseRowStart[jp + 1];
                            m++
                        )
                        {
                            if (ir == coarseColumn[m])
                            {
                                activeCoarseLower[m] += ra*coeffP[indexP];
                                break;
                            }
                        }
                    }
                    else if (ir == jp)
                    {
                        // Found COARSE diagonal
                        activeCoarseDiag[ir] += ra*coeffP[indexP];
                    }
                    else
                    {
                        // Found upper COARSE tringle
                        // Search for coefficient insertion.  HJ: this needs optimisation
                        for
                        (
                            label m = coarseRowStart[ir];
                            m < coarseRowStart[ir + 1];
                            m++
                        )
                        {
                            if (jp == coarseColumn[m])
                            {
                                activeCoarseUpper[m] += ra*coeffP[indexP];
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "autoPtr<BlockAMGLevel<Type> > "
            "BlockMatrixAgglomeration<Type>::restrictMatrix() const"
        )   << "Matrix diagonal of scalar or linear type not implemented"
            << abort(FatalError);
    }
    Info<< "End matrix restriction.  Level size: " << nCoarseEqns_
        << endl;
    // Create and return BlockAMGLevel
    return autoPtr<BlockAMGLevel<Type> >
    (
        new coarseBlockAMGLevel<Type>
        (
            coarseAddrPtr,
            coarseMatrixPtr,
            this->dict(),
            this->type(),
            0,                      // Group size is not used in SAMG
            this->minCoarseEqns()
        )
    );
}


template<class Type>
void Foam::BlockMatrixSelection<Type>::restrictResidual
(
    const Field<Type>& res,
    Field<Type>& coarseRes
) const
{
    coarseRes = pTraits<Type>::zero;

    // Get reference to restriction matrix
    const crMatrix& R = *Rptr_;

    // Get restriction addressing
    const crAddressing& crR = R.crAddr();
    const labelList& rowStartR = crR.row();
    const labelList& colR = crR.col();

    // Coefficients of restriction
    const scalarField& coeffR = R.coeffs();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    // Multiply the residual with restriction weights to obtain the initial
    // coarse residual

    for (label i = 0; i < nCoarseEqns_; i++)
    {
        for (label index = rowStartR[i]; index < rowStartR[i + 1]; index++)
        {
            // Multiply each coeff in row of restriction with the corresponding
            // residual coefficient (column index of R is the same as row index
            // of residual)

            // Column of R
            label jr = colR[index];

            coarseRes[i] += mult(coeffR[index], res[jr]);
        }
    }
}

template<class Type>
void Foam::BlockMatrixSelection<Type>::prolongateCorrection
(
    Field<Type>& x,
    const Field<Type>& coarseX
) const
{
    // Get reference to prolongation matrix
    const crMatrix& P = *Pptr_;

    // Get prolongation addressing
    const crAddressing& crP = P.crAddr();
    const labelList& rowStartP = crP.row();
    const labelList& colP = crP.col();
    const label sizeP = crP.nRows();

    // Coefficients of prolongation
    const scalarField& coeffP = P.coeffs();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    // Multiply the coarse level solution with prolongation and obtain fine
    // level solution

    for (label i = 0; i < sizeP; i++)
    {
        for (label index = rowStartP[i]; index < rowStartP[i + 1]; index++)
        {
            // Multiply each coeff in row of prolongation with the corresponding
            // coarse residual coefficient (column index of prolongation must be
            // the same as row index of coarse residual coeff)

            // Column of P
            label jp = colP[index];

            x[i] += mult(coeffP[index], coarseX[jp]);
        }
    }
}


// ************************************************************************* //
