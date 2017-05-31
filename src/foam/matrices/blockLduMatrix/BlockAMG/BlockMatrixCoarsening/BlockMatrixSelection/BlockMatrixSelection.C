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
    Classical AMG coarsening algorithm for block matrices.

Author
    Tessa Uroic, FMENA

\*----------------------------------------------------------------------------*/

#include "BlockMatrixSelection.H"
#include "coeffFields.H"
#include "BlockSAMGInterfaceField.H"
#include "coarseBlockAMGLevel.H"
#include "PriorityList.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Factor for defining strong negative-coupling of variables
template<class Type>
const Foam::debug::tolerancesSwitch
Foam::BlockMatrixSelection<Type>::epsilon_
(
    "samgStrongConnectionFactor",
    0.25
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::crMatrix>
Foam::BlockMatrixSelection<Type>::filterProlongation
(
    const crMatrix& prolongationMatrix,
    const labelList& fineFaceCells
) const
{
    // Get the addressing and coefficients of the prolongation matrix
    // on my side
    const labelList& pRowStart = prolongationMatrix.crAddr().rowStart();
    const labelList& pColumn = prolongationMatrix.crAddr().column();
    const scalarField& pCoeffs = prolongationMatrix.coeffs();
    const label pNCols = prolongationMatrix.crAddr().nCols();

    // Count how many weighting factors should be sent across - to avoid
    // using dynamic lists
    label countWeights = 0;

    forAll (fineFaceCells, i)
    {
        label end = fineFaceCells[i] + 1;
        label start = fineFaceCells[i];
        countWeights += pRowStart[end] - pRowStart[start];
    }

    // Create filtered prolongation matrix arrays
    labelList filteredRow(fineFaceCells.size() + 1);
    labelList filteredCol(countWeights, 0);
    scalarField filteredCoeffs(countWeights, 0);

    // Select the part of the prolongation matrix to send (only for equations
    // on processor boundary)

    // Mark the start of the compressed row addressing
    // Please note: row start addressing and the owner (faceCells) are linked by
    // their address, i.e. the beginning of row faceCells[i] will be at
    // rowStart[i], and NOT at rowStart[faceCells[i]]!
    filteredRow[0] = 0;

    forAll (fineFaceCells, i)
    {
        label nCoeffs = filteredRow[i];
        const label start = fineFaceCells[i];
        const label end = fineFaceCells[i] + 1;

        // Copy coefficients from the prolongation matrix on my side to the
        // processor boundary (filtered) prolongation matrix
        for (label k = pRowStart[start]; k < pRowStart[end]; k++)
        {
            filteredCoeffs[nCoeffs] = pCoeffs[k];
            filteredCol[nCoeffs] = pColumn[k];
            nCoeffs++;
        }

        // Update row start addressing for next row
        filteredRow[i + 1] = nCoeffs;
    }

    // Create filtered prolongation matrix
    autoPtr<crMatrix> tFilteredP
    (
        new crMatrix
        (
            fineFaceCells.size(),
            pNCols,
            filteredRow,
            filteredCol
        )
    );
    crMatrix& filteredP = tFilteredP();

    // Transfer coefficients
    filteredP.coeffs().transfer(filteredCoeffs);

    return tFilteredP;
}



template<class Type>
void Foam::BlockMatrixSelection<Type>::calcCoarsening()
{
//------------------------------------------------------------------------------
//             MATRIX DATA: ADDRESSING, COEFFICIENTS, COEFF NORMS
//------------------------------------------------------------------------------

    // Get addressing
    const unallocLabelList& rowStart = matrix_.lduAddr().ownerStartAddr();
    const unallocLabelList& losortAddr = matrix_.lduAddr().losortAddr();
    const unallocLabelList& losortStart = matrix_.lduAddr().losortStartAddr();
    const unallocLabelList& column = matrix_.lduAddr().upperAddr();
    const unallocLabelList& row = matrix_.lduAddr().lowerAddr();
    const label matrixSize = matrix_.lduAddr().size();

    // Note: not taking norm magnitudes.  HJ, 28/Feb/2017

//------------------------------------------------------------------------------
//                         CALCULATE COARSENING NORM
//------------------------------------------------------------------------------

    // Calculate norm for diagonal coefficients
    scalarField normDiag(matrixSize);
    coarseningNormPtr_->normalize(normDiag, matrix_.diag());

    // Calculate norm for upper triangle coeffs (magUpper)
    scalarField normUpper(row.size());
    coarseningNormPtr_->normalize(normUpper, matrix_.upper());

    // Calculate norm for lower triangle coeffs (magLower)
    scalarField normLower(column.size());
    coarseningNormPtr_->normalize(normLower, matrix_.lower());

    // Calculate norm magnitudes
    scalarField magNormUpper = mag(normUpper);
    scalarField magNormLower = mag(normLower);

//------------------------------------------------------------------------------
//                        CALCULATE INTERPOLATION NORM
//------------------------------------------------------------------------------

    // Calculate norm for diagonal coefficients (store into magDiag)
    scalarField interpolationNormDiag(matrixSize);
    interpolationNormPtr_->normalize(interpolationNormDiag, matrix_.diag());

    // Calculate norm for upper triangle coeffs (magUpper)
    scalarField interpolationNormUpper(column.size());
    interpolationNormPtr_->normalize(interpolationNormUpper, matrix_.upper());

    // Calculate norm for lower triangle coeffs (magLower)
    scalarField interpolationNormLower(column.size());
    interpolationNormPtr_->normalize(interpolationNormLower, matrix_.lower());

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

    register label j = -1;

    for (register label i = 0; i < matrixSize; i++)
    {
        scalar signDiag = sign(normDiag[i]);

        for (register label k = rowStart[i]; k < rowStart[i + 1]; k++)
        {
            // Do row coefficient
            scalar magAij = mag(min(signDiag*normUpper[k], 0));

            // Upper triangle
            if (magAij > epsLargestNorm[i])
            {
                epsLargestNorm[i] = magAij;
            }

            // Do col coefficient
            scalar magAji = mag(min(signDiag*normLower[k], 0));

            // Lower triangle
            j = column[k];

            if (magAji > epsLargestNorm[j])
            {
                epsLargestNorm[j] = magAji;
            }
        }
    }

    // Multiply largest norm by epsilon.  This is now used below
    epsLargestNorm *= epsilon_();

    // Count strong elements in each row - for rowStart addressing
    // Note: checking magnitude of coeff, which accounts for both strong
    // positive and negative coefficients.  HJ, 28/Feb/2017
    labelList strongCoeffCounter(matrixSize, 0);

    for (register label i = 0; i < matrixSize; i++)
    {
        scalar signDiag = sign(normDiag[i]);

        for (register label k = rowStart[i]; k < rowStart [i + 1]; k++)
        {
            // Do row coefficient
            scalar magAij = mag(min(signDiag*normUpper[k], 0));

            // Upper triangle
            if (magAij > epsLargestNorm[i])
            {
                    strongCoeffCounter[i]++;
            }

            // Do col coefficient
            scalar magAji = mag(min(signDiag*normLower[k], 0));

            // Column of coefficient k in row i
            j = column[k];

            if (magAji > epsLargestNorm[j])
            {
                strongCoeffCounter[j]++;
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
    const labelList& strongRowStart = strong.crAddr().rowStart();
    labelList& strongCol = strong.column();
    scalarField& strongCoeff = strong.coeffs();

    // Counter for counting the strong elements in a row
    labelList counter(matrixSize, 0);

    for (register label i = 0; i < matrixSize; i++)
    {
        scalar signDiag = sign(normDiag[i]);

        for (register label k = rowStart[i]; k < rowStart [i+1]; k++)
        {
            // Check elements in upper triangle
            scalar magAij = mag(min(signDiag*normUpper[k], 0));

            if (magAij > epsLargestNorm[i])
            {
                // Store the strong elements into crMatrix, use counter
                // to count the number of negative strong elements for
                // each row i
                strongCol[strongRowStart[i] + counter[i]] = j;
                strongCoeff[strongRowStart[i] + counter[i]] =
                    interpolationNormUpper[k];

                counter[i]++;
            }

            // Check elements in lower triangle
            scalar magAji = mag(min(signDiag*normLower[k], 0));

            // Column of coefficient k in row i
            j = column[k];

            if (magAji > epsLargestNorm[j])
            {
                strongCol[strongRowStart[j] + counter[j]] = i;
                strongCoeff[strongRowStart[j] + counter[j]] =
                    interpolationNormLower[k];

                counter[j]++;
            }
        }
    }

//------------------------------------------------------------------------------
//           COARSENING: SORT EQUATIONS INTO COARSE AND FINE SUBSETS
//------------------------------------------------------------------------------

    // Transpose the compressed row matrix to use for coarsening
    crAddressing Taddr = strong.crAddr().T();
    const labelList& tRow = Taddr.rowStart();
    const labelList& tCol = Taddr.column();

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
        equationWeight.set(i, strongRowStart[i + 1] - strongRowStart[i]);

        // Label rows without connections as FINE
        if (strongRowStart[i + 1] == strongRowStart[i])
        {
            rowLabel_[i] = FINE;
        }
    }

    // Start counting coarse equations
    nCoarseEqns_ = 0;

    while (!equationWeight.empty())
    {
        // removeHead = return index with largest weight and remove
        label topElement = equationWeight.removeHead();

        if (rowLabel_[topElement] == UNDECIDED)
        {
            rowLabel_[topElement] = nCoarseEqns_;
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

    // From processor boundary
    const typename BlockLduInterfaceFieldPtrsList<Type>::Type& interfaceFields =
        matrix_.interfaces();

    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            // Get norm of boundary coefficients
            scalarField normBou(matrix_.coupleUpper()[intI].size());
            interpolationNormPtr_->normalize
            (
                normBou,
                matrix_.coupleUpper()[intI]
            );

            // Get addressing
            const labelList& owner =
                interfaceFields[intI].coupledInterface().faceCells();

            forAll (normBou, coeffI)
            {
                // Get row from faceCells
                const label i = owner[coeffI];

                // Get sign of diagonal coeff
                scalar signDiag = sign(interpolationNormDiag[i]);

                // Adjust sign of off-diag coeff
                // Note: additional minus because the sign of interface
                // coeffs is opposite from the normal matrix off-diagonal

                scalar coeff = -signDiag*normBou[coeffI];

                // Add negative/positive contribution into corresponding field
                // which contributes to prolongation weight factor
                positiveElemSum += Foam::max(coeff, 0);
                negativeElemSum += Foam::min(coeff,0);
            }
        }
    }


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
            scalar signDiag = sign(interpolationNormDiag[i]);

            // Upper triangle - use rowStart
            for (label k = rowStart[i]; k < rowStart[i + 1]; k++)
            {
                scalar coeff = signDiag*interpolationNormUpper[k];

                positiveElemSum[i] += Foam::max(coeff, 0);
                negativeElemSum[i] += Foam::min(coeff, 0);
            }

            // Lower triangle - use losort to go row-wise (Why? This
            // is a pain)
            // Because I have to compare the sign to the
            // corresponding diagonal element in row i.
            for (label k = losortStart[i]; k < losortStart[i + 1]; k++)
            {
                label index = losortAddr[k];

                scalar coeff = signDiag*interpolationNormLower[index];

                positiveElemSum[i] += Foam::max(coeff, 0);
                negativeElemSum[i] += Foam::min(coeff, 0);
            }
        }
    }

    // Adding only strong coarse contributions for each row (denominator
    // of scaling factor)

    for (label i = 0; i < matrixSize; i++)
    {
        if (rowLabel_[i] == FINE)
        {
            scalar signDiag = sign(interpolationNormDiag[i]);

            for
            (
                label index = strongRowStart[i];
                index < strongRowStart[i + 1];
                index++
            )
            {
                // Get column of strong coeff
                label j = strongCol[index];
                scalar coeff = signDiag*strongCoeff[index];

                if (rowLabel_[j] != FINE)
                {
                    posCoarseElemSum[i] += Foam::max(coeff, 0);
                    negCoarseElemSum[i] += Foam::min(coeff, 0);
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
            scalar diagonalCoeff = interpolationNormDiag[i];

            // Check whether this Eqn has COARSE negative strong
            // contributions!

            if (negCoarseElemSum[i] == 0)
            {
                //#   ifdef FULLDEBUG

                // Positive connections are used for interpolation! Diagonal
                // dominance is not conserved!
                Info << "Interpolation from positive neighbours!" << nl
                     << "Equation " << i << endl;

                //#   endif

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

                    if (rowLabel_[j] != FINE)
                    {
                        prolongationCoeff.append
                        (
                            (-positiveElemSum[i])/posCoarseElemSum[i]*
                            strongCoeff[k]/(diagonalCoeff)
                        );

                        prolongationCol.append(rowLabel_[j]);
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
                        rowLabel_[j] != FINE &&
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

                        prolongationCol.append(rowLabel_[j]);
                        rowCount++;
                    }
                }

            }
        }
        else // The row is COARSE - injection
        {
            prolongationCoeff.append(1);
            prolongationCol.append(rowLabel_[i]);
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

    if (blockLduMatrix::debug >= 2)
    {
        Pout<< "Coarse level size: " << nCoarseEqns_;

        if (coarsen_)
        {
            Pout<< ".  Accepted" << endl;
        }
        else
        {
            Pout<< ".  Rejected" << endl;
        }
    }

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
    coarseningNormPtr_
    (
        BlockCoeffNorm<Type>::New(dict.subDict("coarseningNorm"))
    ),
    interpolationNormPtr_
    (
        BlockCoeffNorm<Type>::New(dict.subDict("interpolationNorm"))
    ),
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
    const labelList& rowStartR = crR.rowStart();
    const labelList& colR = crR.column();

    // Matrix A addressing
    const unallocLabelList& rowStartA = matrix_.lduAddr().ownerStartAddr();
    const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();
    const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
    const unallocLabelList& losortAddr = matrix_.lduAddr().losortAddr();
    const unallocLabelList& losortStart = matrix_.lduAddr().losortStartAddr();

    // Prolongation addressing
    const labelList& rowStartP = crP.rowStart();
    const labelList& colP = crP.column();

    // coeffLabel is used for checking if a contribution in that address already
    // exists
    labelList coeffLabel(nCoarseEqns_, -1);
    label nCoarseCoeffs = 0;

    // Create coarse addressing - owner and neighbour pairs
    // HJ: give a better guess of size of coarse off-diagonal
    HashTable<label, labelPair, Hash<labelPair> > coarseOwnNbr;
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
        coeffLabel = -1;

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

                            coarseOwnNbr.insert
                            (
                                labelPair(ir, jp),
                                nCoarseCoeffs
                            );

                            coarseOwner.append(ir);
                            coarseNeighbour.append(jp);

                            nCoarseCoeffs++;
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

                            coarseOwnNbr.insert
                            (
                                labelPair(ir, jp),
                                nCoarseCoeffs
                            );

                            coarseOwner.append(ir);
                            coarseNeighbour.append(jp);

                            nCoarseCoeffs++;
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

                        coarseOwnNbr.insert
                        (
                            labelPair(ir, jp),
                            nCoarseCoeffs
                        );

                        coarseOwner.append(ir);
                        coarseNeighbour.append(jp);

                        nCoarseCoeffs++;
                    }
                }
            }
        }
    }

//------------------------------------------------------------------------------
//                      CREATE COARSE MATRIX ADDRESSING
//------------------------------------------------------------------------------

    // Set size of dynamic list
    coarseOwner.setSize(nCoarseCoeffs);
    coarseNeighbour.setSize(nCoarseCoeffs);

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

    const typename BlockLduInterfaceFieldPtrsList<Type>::Type& interfaceFields =
        matrix_.interfaces();

    // Set the coarse interfaces and coefficients
    lduInterfacePtrsList coarseInterfaces(interfaceFields.size());

    // Initialise transfer of matrix prolongation on the interface

    // Store owner and neighbour prolongation to avoid tangled communications
    PtrList<crMatrix> ownInterfaceProlongation(interfaceFields.size());
    PtrList<crMatrix> nbrInterfaceProlongation(interfaceFields.size());

    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            const labelList& fineFaceCells =
                interfaceFields[intI].coupledInterface().faceCells();

            // Filter local prolongation matrix and return to
            // ownInterfaceProlongation
            ownInterfaceProlongation.set
            (
                intI,
                filterProlongation(P, fineFaceCells).ptr()
            );

            interfaceFields[intI].coupledInterface().initProlongationTransfer
            (
                Pstream::blocking,
                ownInterfaceProlongation[intI]
            );
        }
    }

    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            const lduInterface& fineInterface =
                interfaceFields[intI].coupledInterface();

            nbrInterfaceProlongation.set
            (
                intI,
                fineInterface.prolongationTransfer
                (
                    Pstream::blocking,
                    ownInterfaceProlongation[intI]
                ).ptr()
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
                SAMGInterface::New
                (
                    coarseAddrPtr(),
                    P,
                    coarseInterfaces,
                    fineInterface,
                    nbrInterfaceProlongation[intI]
                ).ptr()
            );
        }
    }

    labelListList coarseInterfaceAddr(interfaceFields.size());

    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            const SAMGInterface& coarseInterface =
                refCast<const SAMGInterface>(coarseInterfaces[intI]);

            coarseInterfaceAddr[intI] = coarseInterface.faceCells();
        }
    }

    // Add interfaces
    coarseAddrPtr->addInterfaces
    (
        coarseInterfaces,
        coarseInterfaceAddr,
        matrix_.patchSchedule()
    );


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
                            label face = coarseOwnNbr[labelPair(jp, ir)];
                            activeCoarseLower[face] += ra*coeffP[indexP];
                        }
                        else if (ir == jp)
                        {
                            // Found COARSE diagonal coefficient
                            activeCoarseDiag[ir] += ra*coeffP[indexP];
                        }
                        else
                        {
                            // Found upper COARSE triangle
                            label face = coarseOwnNbr[labelPair(ir, jp)];
                            activeCoarseUpper[face] += ra*coeffP[indexP];
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
                            label face = coarseOwnNbr[labelPair(jp, ir)];
                            activeCoarseLower[face] += ra*coeffP[indexP];
                        }
                        else if (ir == jp)
                        {
                            // Found COARSE diagonal
                            activeCoarseDiag[ir] += ra*coeffP[indexP];
                        }
                        else
                        {
                            // Found upper COARSE triangle
                            label face = coarseOwnNbr[labelPair(ir, jp)];
                            activeCoarseUpper[face] += ra*coeffP[indexP];
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
                        label face = coarseOwnNbr[labelPair(jp, ir)];
                        activeCoarseLower[face] += ra*coeffP[indexP];
                    }
                    else if (ir == jp)
                    {
                        // Found COARSE diagonal
                        activeCoarseDiag[ir] += ra*coeffP[indexP];
                    }
                    else
                    {
                        // Found upper COARSE triangle
                        label face = coarseOwnNbr[labelPair(ir, jp)];
                        activeCoarseUpper[face] += ra*coeffP[indexP];
                    }
                }
            }
        }

        // Get interfaces from coarse matrix
        typename BlockLduInterfaceFieldPtrsList<Type>::Type&
            coarseInterfaceFieldsTransfer = coarseMatrix.interfaces();

        // Aggolmerate the upper and lower coupled coefficients
        forAll (interfaceFields, intI)
        {
            if (interfaceFields.set(intI))
            {
                const SAMGInterface& coarseInterface =
                    refCast<const SAMGInterface>(coarseInterfaces[intI]);

                coarseInterfaceFieldsTransfer.set
                (
                    intI,
                    BlockSAMGInterfaceField<Type>::New
                    (
                        coarseInterface,
                        interfaceFields[intI]
                    ).ptr()
                );

                // Since the type of clustering is now templated, clustering
                // of block coefficients must be done by a FIELD (not interface)
                // via a new set of virtual functions
                // HJ, 16/Mar/2016

                // Note: in the scalar AMG, clustering is done by the interface
                // (always scalar) but in the block matrix it is done by a
                // templated block interface field
                // HJ, 16/Mar/2016

                // Cast the interface into AMG type
                const BlockSAMGInterfaceField<Type>& coarseField =
                    refCast<const BlockSAMGInterfaceField<Type> >
                    (
                        coarseInterfaceFieldsTransfer[intI]
                    );

                coarseMatrix.coupleUpper().set
                (
                    intI,
                    coarseField.selectBlockCoeffs
                    (
                        matrix_.coupleUpper()[intI]
                    )
                );

                coarseMatrix.coupleLower().set
                (
                    intI,
                    coarseField.selectBlockCoeffs
                    (
                        matrix_.coupleLower()[intI]
                    )
                );
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
    const labelList& rowStartR = crR.rowStart();
    const labelList& colR = crR.column();

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
    const labelList& rowStartP = crP.rowStart();
    const labelList& colP = crP.column();
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
