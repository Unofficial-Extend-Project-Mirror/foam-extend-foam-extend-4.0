/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

// Factor defining strong negative-coupling of variables
template<class Type>
const Foam::debug::tolerancesSwitch
Foam::BlockMatrixSelection<Type>::epsilon_
(
    "blockSamgEpsilon",
    0.2
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
    const label nRows = matrix_.lduAddr().size();
    const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();
    const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
    const unallocLabelList& losortAddr = matrix_.lduAddr().losortAddr();
    const unallocLabelList& ownerStart = matrix_.lduAddr().ownerStartAddr();
    const unallocLabelList& losortStart = matrix_.lduAddr().losortStartAddr();

    // Note: not taking norm magnitudes.  HJ, 28/Feb/2017

//------------------------------------------------------------------------------
//                              CALCULATE NORM
//------------------------------------------------------------------------------

    // Calculate norm for diagonal coefficients
    scalarField normDiag(nRows);
    normPtr_->normalize(normDiag, matrix_.diag());

    // Calculate norm for upper triangle coeffs (magUpper)
    scalarField normUpper(upperAddr.size());
    normPtr_->normalize(normUpper, matrix_.upper());

    // Calculate norm for lower triangle coeffs (magLower)
    scalarField normLower(upperAddr.size());
    normPtr_->normalize(normLower, matrix_.lower());

    // Calculate norm magnitudes
    scalarField magNormDiag = mag(normDiag);
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
    scalarField epsilonStrongCoeff(nRows, 0);

    // Select the strongest coefficient in each row
    for (label i = 0; i < nRows; i++)
    {
        const scalar signDiag = sign(normDiag[i]);

        // Do lower triangle coefficient for the row first
        for (label jp = losortStart[i]; jp < losortStart[i + 1]; jp++)
        {
            const scalar magAij =
                mag(min(signDiag*normLower[losortAddr[jp]], 0));

            if (magAij > epsilonStrongCoeff[i])
            {
                epsilonStrongCoeff[i] = magAij;
            }
        }

        // Do upper triangle coefficient for the row
        for (label ip = ownerStart[i]; ip < ownerStart[i + 1]; ip++)
        {
            const scalar magAij = mag(min(signDiag*normUpper[ip], 0));

            if (magAij > epsilonStrongCoeff[i])
            {
                epsilonStrongCoeff[i] = magAij;
            }
        }
    }

    // Multiply largest norm by epsilon.  This is now used below
    epsilonStrongCoeff *= epsilon_();

    // Count strong elements in each row - for row addressing
    // Note: checking magnitude of coeff, which accounts for both strong
    // positive and negative coefficients.  HJ, 28/Feb/2017
    labelList strongCoeffCounter(nRows, 0);

    for (label i = 0; i < nRows; i++)
    {
        const scalar signDiag = sign(normDiag[i]);

        // Do lower triangle coefficient for the row first
        for (label jp = losortStart[i]; jp < losortStart[i + 1]; jp++)
        {
            const scalar magAij =
                mag(min(signDiag*normLower[losortAddr[jp]], 0));

            if (magAij > epsilonStrongCoeff[i])
            {
                strongCoeffCounter[i]++;
            }
        }

        // Do upper triangle coefficient for the row
        for (label ip = ownerStart[i]; ip < ownerStart[i + 1]; ip++)
        {
            const scalar magAij = mag(min(signDiag*normUpper[ip], 0));

            if (magAij > epsilonStrongCoeff[i])
            {
                strongCoeffCounter[i]++;
            }
        }
    }

    // Create a crMatrix that will store all of the strong elements of each row
    // (some will become FINE and some COARSE - and this will determine the
    // direct and standard interpolation procedures)

    crMatrix strong(nRows, nRows, strongCoeffCounter);

    // Set addressing for the matrix:
    // stongCol and strongElement are arrays needed to create a compressed row
    // matrix
    const labelList& strongRow = strong.crAddr().rowStart();
    labelList& strongCol = strong.column();
    scalarField& strongCoeff = strong.coeffs();

    // Reset strongCoeffCounter for re-use
    strongCoeffCounter = 0;

    for (label i = 0; i < nRows; i++)
    {
        const scalar signDiag = sign(normDiag[i]);

        // Do lower triangle coefficient for the row first
        for (label jp = losortStart[i]; jp < losortStart[i + 1]; jp++)
        {
            const scalar magAij =
                mag(min(signDiag*normLower[losortAddr[jp]], 0));

            if (magAij > epsilonStrongCoeff[i])
            {
                strongCol[strongRow[i] + strongCoeffCounter[i]] =
                    lowerAddr[losortAddr[jp]];
                strongCoeff[strongRow[i] + strongCoeffCounter[i]] =
                    normLower[losortAddr[jp]];

                strongCoeffCounter[i]++;
            }
        }

        // Do upper triangle coefficient for the row
        for (label ip = ownerStart[i]; ip < ownerStart[i + 1]; ip++)
        {
            const scalar magAij = mag(min(signDiag*normUpper[ip], 0));

            if (magAij > epsilonStrongCoeff[i])
            {
                strongCol[strongRow[i] + strongCoeffCounter[i]] = upperAddr[ip];
                strongCoeff[strongRow[i] + strongCoeffCounter[i]] =
                    normUpper[ip];

                strongCoeffCounter[i]++;
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
    rowLabel_.setSize(nRows, UNDECIDED);

    PriorityList<label> equationWeight(nRows);

    for (label i = 0; i < nRows; i++)
    {
        // Set weights for each equation
        // Count equations that my equation is strongly influencing
        // (dependants). TU, 7/Jul/2017
        equationWeight.set(i, tRow[i + 1] - tRow[i]);
    }

    // Mark disconnected rows as fine. This also removes solo
    // equations from coarsening
    // HJ, 29/Jul/2017
    for (label i = 0; i < nRows; i++)
    {
        // Label rows without strong connections as FINE
        if (strongRow[i + 1] == strongRow[i])
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
            // Make highest equation choice coarse
            rowLabel_[topElement] = nCoarseEqns_;
            nCoarseEqns_++;

            // Decrement weights of neighbouring equations
            for
            (
                label ip = strongRow[topElement];
                ip < strongRow[topElement + 1];
                ip++
            )
            {
                label j = strongCol[ip];

                if (rowLabel_[j] == UNDECIDED)
                {
                    equationWeight.updateWeight
                    (
                        j,
                        equationWeight.weights()[j] - 1
                    );
                }
            }

            // Make all neighbours fine and increment weight
            for
            (
                label ip = tRow[topElement];
                ip < tRow[topElement + 1];
                ip++
            )
            {
                label j = tCol[ip];

                if (rowLabel_[j] == UNDECIDED)
                {
                    rowLabel_[j] = FINE;

                    for (label jp = strongRow[j]; jp < strongRow[j + 1]; jp++)
                    {
                        label kp = strongCol[jp];

                        if (rowLabel_[kp] == UNDECIDED)
                        {
                            equationWeight.updateWeight
                            (
                                kp,
                                equationWeight.weights()[kp] + 2
                            );
                        }
                    }
                }
            }
        }
    }

//------------------------------------------------------------------------------
//              CALCULATING CONTRIBUTIONS TO THE SCALING FACTOR
//------------------------------------------------------------------------------

    // Sum of negative elements in row
    scalarField num(nRows, 0);

    // Sum of positive elements in row (sign equal to diagonal)
    scalarField Dii(sign(normDiag)*normDiag);

    scalar Dij, den, signDii;

    // Coupled boundary tratment
    // In preparation for the calculation of interpolation matrices
    // take out all coupled processor boundary coeffs from the
    // equations near the boundary
    // Note:
    // bouCoeffs() = upper triangular coupled coeff
    // intCoeffs() = lower triangular coupled coeff
    // Note:
    // The sign of bouCoeffs and intCoeffs is opposite
    // of the sign of internal off-diag coeffs (they are on rhs)
    // Note:
    // On the local processor, we only account for upper coeffs, as the
    // column coeff is handled by the other processor
    // HJ, 16/May/2017

    // Naming convention for off-diagonal coefficients
    // Positive coefficient = same sign as diagonal (bad)
    // Negative coefficient = opposite sign od diagonal (good)


    // Collect contributions from coupled boundaries
    const typename BlockLduInterfaceFieldPtrsList<Type>::Type& interfaceFields =
        matrix_.interfaces();

    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            // Get norm of boundary coefficients
            scalarField normCplUpper(matrix_.coupleUpper()[intI].size());
            normPtr_->normalize
            (
                normCplUpper,
                matrix_.coupleUpper()[intI]
            );

            // Get addressing
            const labelList& faceCells =
                interfaceFields[intI].coupledInterface().faceCells();

            forAll (normCplUpper, coeffI)
            {
                // Get row from faceCells
                const label i = faceCells[coeffI];

                // Get sign of diagonal coeff
                scalar signDiag = sign(normDiag[i]);

                // Adjust sign of off-diag coeff
                // Note: additional minus because the sign of interface
                // coeffs is opposite from the normal matrix off-diagonal
                Dij = -signDiag*normCplUpper[coeffI];

                // Add negative/positive contribution into corresponding field
                // which contributes to prolongation weight factor
                Dii[i] += Foam::max(Dij, 0);
                num[i] += Foam::min(Dij, 0);   // HJ, HERE!!!
            }
        }
    }


//------------------------------------------------------------------------------
//                    CALCULATING PROLONGATION WEIGHTS
//------------------------------------------------------------------------------

    // Cannot create P addressing in first pass.  Count and resize arrays
    label maxPCount = 0;

    forAll (rowLabel_, i)
    {
        if (rowLabel_[i] == FINE)
        {
            // Interpolation involves neighbourhood
            maxPCount += strongRow[i + 1] - strongRow[i];
        }
        else
        {
            // Coarse point, injection
            maxPCount++;
        }
    }

    labelList pRow(nRows + 1);
    labelList pCol(maxPCount, 0);
    scalarField pCoeff(maxPCount, 0);

    // Start row assembly
    pRow[0] = 0;

    for (label i = 0; i < nRows; i++)
    {
        label rowCount = pRow[i];

        // Handle the unknown sign of diagonal: multiplying
        // row coeffs by the sign
        signDii = sign(normDiag[i]);

        // Do lower triangle coefficient for the row first
        for (label jp = losortStart[i]; jp < losortStart[i + 1]; jp++)
        {
            // Adjust sign of off-diag coeff
            Dij = sign(normDiag[i])*normLower[losortAddr[jp]];

            // Add negative coeff to num
            // Add positive coeff to diag to eliminate it
            num[i] += Foam::min(Dij, 0);
            Dii[i] += Foam::max(Dij, 0);
        }

        // Do upper triangle coefficient for the row
        for (label ip = ownerStart[i]; ip < ownerStart[i + 1]; ip++)
        {
            // Adjust sign of off-diag coeff
            Dij = signDii*normUpper[ip];

            // Add negative coeff to num
            // Add positive coeff to diag to eliminate it
            num[i] += Foam::min(Dij, 0);
            Dii[i] += Foam::max(Dij, 0);
        }

        // Row i completed.  Calculate weights
        if (rowLabel_[i] == FINE)
        {
            // Fine equation
            den = 0;

            for
            (
                label sip = strongRow[i];
                sip < strongRow[i + 1];
                sip++
            )
            {
                const label js = strongCol[sip];

                if (rowLabel_[js] != FINE)
                {
                    den += strongCoeff[sip];
                }
            }

            for
            (
                label sip = strongRow[i];
                sip < strongRow[i + 1];
                sip++
            )
            {
                const label js = strongCol[sip];

                if (rowLabel_[js] != FINE)
                {
                    // PROLONGATION NOT SUMMING INTO 1 FOR DIAGONALLY UNEQUAL
                    // ROWS:
                    // Prolongation coefficient with scaling
                    pCoeff[rowCount] = -(num[i]/den)*strongCoeff[sip]/Dii[i];

                    pCol[rowCount] = rowLabel_[js];
                    rowCount++;
                }
            }
        }
        else
        {
            // Coarse equation
            pCoeff[rowCount] = 1;
            pCol[rowCount] = rowLabel_[i];
            rowCount++;
        }

        // Grab row start/end
        pRow[i + 1] = rowCount;
    }

    // Resize column and coeffs
    pCoeff.setSize(pRow[nRows]);
    pCol.setSize(pRow[nRows]);

    // Create prolongation matrix
    Pptr_ = new crMatrix(nRows, nCoarseEqns_, pRow, pCol);
    crMatrix& prolongation = *Pptr_;

    prolongation.coeffs().transfer(pCoeff);

    // Check prolongation matrix
    if (blockLduMatrix::debug > 2)
    {
        scalarField sumRow(nRows, 0);
        const labelList& prolongationRow = Pptr_->crAddr().rowStart();
        const scalarField& prolongationCoeff = Pptr_->coeffs();

        for (label rowI = 0; rowI < nRows; rowI++)
        {
            for
            (
                label colI = prolongationRow[rowI];
                colI < prolongationRow[rowI + 1];
                colI++
            )
            {
                sumRow[rowI] += prolongationCoeff[colI];
            }
        }

        scalar minSumRow = min(sumRow);

        if (minSumRow < 0.99)
        {
            Pout<< "sumRow (min, max) = (" << min(sumRow) << " "
                << max(sumRow) << ")" << endl;
        }
    }

    // The decision on parallel agglomeration needs to be made for the
    // whole gang of processes; otherwise I may end up with a different
    // number of agglomeration levels on different processors.
    if
    (
        nCoarseEqns_ > this->minCoarseEqns()
     && 3*nCoarseEqns_ <= 2*nRows
    )
    {
        coarsen_ = true;
    }

    reduce(coarsen_, andOp<bool>());

    if (blockLduMatrix::debug >= 3)
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
    normPtr_(BlockCoeffNorm<Type>::New(dict)),
    nCoarseEqns_(0),
    coarsen_(false),
    Pptr_(nullptr),
    Rptr_(nullptr),
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
    const labelList& rowR = crR.rowStart();
    const labelList& colR = crR.column();

    // Matrix A addressing
    const unallocLabelList& rowA = matrix_.lduAddr().ownerStartAddr();
    const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();

    // Addressing for lower triangle loop
    const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();
    const unallocLabelList& losortAddr = matrix_.lduAddr().losortAddr();
    const unallocLabelList& losortStart = matrix_.lduAddr().losortStartAddr();

    // Prolongation addressing
    const labelList& rowP = crP.rowStart();
    const labelList& colP = crP.column();

    // In order to avoid searching for the off-diagonal coefficient,
    // a mark array is used for each row's assembly.
    // Mark records the index of the off-diagonal coefficient in the row
    // for each neighbour entry.  It is reset after completing each row of the
    // coarse matrix.
    // HJ, 28/Apr/2017
    labelList coeffLabel(nCoarseEqns_, -1);

    // Create coarse addressing: record neighbours for every row
    // Note: neighbourt indices are added as they appear and will be sorted
    // on completion
    List<labelHashSet> coarseNbrsSets(nCoarseEqns_);

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
    //         Multiply matrix A coeff with R coeff if rowA == colR
    //         The address of the resulting coeff is (rowR, colA)
    //     RA*P multiplication
    //         Multiply P coeff with RA coeff if rowP == colRA
    //         The address of the resulting coeff is (rowRA, colP), that is
    //         the address in the COARSE matrix is (rowR, colP)

    // NOTE: Letters i and j are used to denote the row and the col of the
    // matrix, respectively. In addition, row of matrix A is denoted ia, col
    // ja, row of prolongation is ip, etc.

    // Loop through rows of R
    for (label ir = 0; ir < nCoarseEqns_; ir++)
    {
        // Compressed row format - get indices of coeffsR in row ir
        for (label indexR = rowR[ir]; indexR < rowR[ir + 1]; indexR++)
        {
            // Col of R coeff
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

                // Col of coeff in A
                const label ja = lowerAddr[faceA];

                // Go into the corresponding row of prolongation to find
                // contributions
                for
                (
                    label indexP = rowP[ja];
                    indexP < rowP[ja + 1];
                    indexP++
                )
                {
                    // Col of coeff in P
                    const label jp = colP[indexP];

                    // Check the address of the coefficient in COARSE matrix
                    // (ir = row, jp = col)
                    if (ir > jp)
                    {
                        // This coeff belongs to the lower triangle
                        // Record it in the upper triangle
                        coarseNbrsSets[jp].insert(ir);
                    }
                    else if (ir == jp)
                    {
                        // Diag coeff ignored
                    }
                    else
                    {
                        // This coeff belongs to the upper triangle
                        // Record it in the upper triangle
                        coarseNbrsSets[ir].insert(jp);
                    }
                }
            }

            // UPPER TRIANGLE
            for
            (
                label faceA  = rowA[jr];
                faceA < rowA[jr + 1];
                faceA++
            )
            {
                // Get col of coeff in A
                const label ja = upperAddr[faceA];

                // Go into the corresponding row of prolongation to find
                // contributions
                for
                (
                    label indexP = rowP[ja];
                    indexP < rowP[ja + 1];
                    indexP++
                )
                {
                    // Get col of coefficient in P
                    const label jp = colP[indexP];

                    // Check the address of the coefficient in COARSE matrix
                    // (ir = row, jp = col)
                    if (ir > jp)
                    {
                        // This coeff belongs to the lower triangle
                        // Record it in the upper triangle
                        coarseNbrsSets[jp].insert(ir);
                    }
                    else if (ir == jp)
                    {
                        // Diag coeff ignored
                    }
                    else
                    {
                        // This coeff belongs to the upper triangle
                        // Record it in the upper triangle
                        coarseNbrsSets[ir].insert(jp);
                    }
                }
            }

            // DIAGONAL
            for
            (
                label indexP = rowP[jr];
                indexP < rowP[jr + 1];
                indexP++
            )
            {
                // Col of coefficient in P
                const label jp = colP[indexP];

                // Check the address of the coefficient
                // (ir = row, jp = col)
                if (ir > jp)
                {
                    // This coeff belongs to the lower triangle
                    // Record it in the upper triangle
                    coarseNbrsSets[jp].insert(ir);
                }
                else if (ir == jp)
                {
                    // Diag coeff ignored
                }
                else
                {
                    // This coeff belongs to the upper triangle
                    // Record it in the upper triangle
                    coarseNbrsSets[ir].insert(jp);
                }
            }
        }
    }

//------------------------------------------------------------------------------
//                      CREATE COARSE MATRIX ADDRESSING
//------------------------------------------------------------------------------

    // Count coarse coeffs
    label nCoarseCoeffs = 0;
    forAll (coarseNbrsSets, rowI)
    {
        nCoarseCoeffs += coarseNbrsSets[rowI].size();
    }

    // Create owner and neighbour lists
    labelList coarseOwner(nCoarseCoeffs);
    labelList coarseNeighbour(nCoarseCoeffs);

    // Fill the owner and neighbour lists from sets, with sorting
    label coeffI = 0;

    forAll (coarseNbrsSets, rowI)
    {
        const labelList curNbrs = coarseNbrsSets[rowI].sortedToc();

        forAll (curNbrs, nbrI)
        {
            coarseOwner[coeffI] = rowI;
            coarseNeighbour[coeffI] = curNbrs[nbrI];
            coeffI++;
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
            false           // DO NOT REUSE STORAGE!
        )
    );

//------------------------------------------------------------------------------
//                       CREATE COARSE MATRIX INTERFACES
//------------------------------------------------------------------------------

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

            // Use filtered prolongation on master and slave side
            // HJ, 7/Jun/2017
            coarseInterfaces.set
            (
                intI,
                SAMGInterface::New
                (
                    coarseAddrPtr(),
                    ownInterfaceProlongation[intI],
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

    // Add interfaces to coarse matrix addressing
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

    // Build matrix coefficients
    this->updateMatrix(coarseMatrix);

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
void Foam::BlockMatrixSelection<Type>::updateMatrix
(
    BlockLduMatrix<Type>& coarseMatrix
) const
{
    // Get references to restriction and prolongation matrix
    const crMatrix& R = *Rptr_;
    const crMatrix& P = *Pptr_;

    // Get connectivity
    const crAddressing& crR = R.crAddr();
    const crAddressing& crP = P.crAddr();

    // Restriction addressing
    const labelList& rowR = crR.rowStart();
    const labelList& colR = crR.column();

    // Matrix A addressing
    const unallocLabelList& rowA = matrix_.lduAddr().ownerStartAddr();
    const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();

    // Get interfaces from fine matrix
    const typename BlockLduInterfaceFieldPtrsList<Type>::Type&
        interfaceFields = matrix_.interfaces();

    // Addressing for lower triangle loop
    const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();
    const unallocLabelList& losortAddr = matrix_.lduAddr().losortAddr();
    const unallocLabelList& losortStart = matrix_.lduAddr().losortStartAddr();

    // Prolongation addressing
    const labelList& rowP = crP.rowStart();
    const labelList& colP = crP.column();

    // In order to avoid searching for the off-diagonal coefficient,
    // a mark array is used for each row's assembly.
    // Mark records the index of the off-diagonal coefficient in the row
    // for each neighbour entry.  It is reset after completing each row of the
    // coarse matrix.
    // HJ, 28/Apr/2017
    labelList coeffLabel(nCoarseEqns_, -1);

    typedef CoeffField<Type> TypeCoeffField;

    TypeCoeffField& coarseUpper = coarseMatrix.upper();
    TypeCoeffField& coarseDiag = coarseMatrix.diag();
    TypeCoeffField& coarseLower = coarseMatrix.lower();

    // Get the coarse interfaces and coefficients
    lduInterfacePtrsList coarseInterfaces = coarseMatrix.mesh().interfaces();

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

    // NOTE: Letters i and j are used to denote the row and the col of the
    // matrix, respectively. In addition, row of matrix A is denoted ia, col
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

        // Get coarse matrix addressing
        const unallocLabelList& rowC = coarseMatrix.lduAddr().ownerStartAddr();
        const unallocLabelList& upperCoarseAddr =
            coarseMatrix.lduAddr().upperAddr();
        const unallocLabelList& lowerCoarseAddr =
            coarseMatrix.lduAddr().lowerAddr();
        const unallocLabelList& losortCoarseAddr =
            coarseMatrix.lduAddr().losortAddr();
        const unallocLabelList& losortCoarseStart =
            coarseMatrix.lduAddr().losortStartAddr();

        // Re-initialise coeffLabel vector
        coeffLabel = -1;

        // Loop through rows of R
        for (label ir = 0; ir < nCoarseEqns_; ir++)
        {
            // Doing new row: expand nbr cells into coeffLabel

            // Upper triangle expand
            for (label coarseK = rowC[ir]; coarseK < rowC[ir + 1]; coarseK++)
            {
                coeffLabel[upperCoarseAddr[coarseK]] = coarseK;
            }

            // Lower triangle expand
            for
            (
                label indexC = losortCoarseStart[ir];
                indexC < losortCoarseStart[ir + 1];
                indexC++
            )
            {
                coeffLabel[lowerCoarseAddr[losortCoarseAddr[indexC]]] =
                    losortCoarseAddr[indexC];
            }

            // Compressed row format, get indices of coeffsR in row ir
            for (label indexR = rowR[ir]; indexR < rowR[ir + 1]; indexR++)
            {
                // Col of coeff in R
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

                    // Col index of coeff in A
                    const label ja = lowerAddr[faceA];

                    // Multiply coefficients of R and A
                    // Address is (rowR, colA)
                    const squareType ra = activeLower[faceA]*coeffR[indexR];

                    // Go into corresponding row of prolongation to find
                    // contributions
                    for
                    (
                        label indexP = rowP[ja];
                        indexP < rowP[ja + 1];
                        indexP++
                    )
                    {
                        // Col of coeff in P
                        const label jp = colP[indexP];

                        // Check the address of the coefficient
                        // (ir = row, jp = col)
                        if (ir > jp)
                        {
                            // Found lower COARSE triangle
                            label face = coeffLabel[jp];
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
                            label face = coeffLabel[jp];
                            activeCoarseUpper[face] += ra*coeffP[indexP];
                        }
                    }
                }

                // UPPER TRIANGLE
                for
                (
                    label faceA = rowA[jr];
                    faceA < rowA[jr + 1];
                    faceA++
                )
                {
                    // Get col of coeff in A
                    const label ja = upperAddr[faceA];

                    const squareType ra = coeffR[indexR]*activeUpper[faceA];

                    for
                    (
                        label indexP = rowP[ja];
                        indexP < rowP[ja + 1];
                        indexP++
                    )
                    {
                        // Get col of coefficient in P
                        const label jp = colP[indexP];

                        if (ir > jp)
                        {
                            // Found lower COARSE triangle
                            label face = coeffLabel[jp];
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
                            label face = coeffLabel[jp];
                            activeCoarseUpper[face] += ra*coeffP[indexP];
                        }
                    }
                }

                // DIAGONAL
                const squareType ra = coeffR[indexR]*activeDiag[jr];

                for
                (
                    label indexP = rowP[jr];
                    indexP < rowP[jr + 1];
                    indexP++
                )
                {
                    const label jp = colP[indexP];

                    if (ir > jp)
                    {
                        // Found lower COARSE triangle
                        label face = coeffLabel[jp];
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
                        label face = coeffLabel[jp];
                        activeCoarseUpper[face] += ra*coeffP[indexP];
                    }
                }
            }

            // Clean out coeffLabel
            // Upper triangle expand
            for (label coarseK = rowC[ir]; coarseK < rowC[ir + 1]; coarseK++)
            {
                coeffLabel[upperCoarseAddr[coarseK]] = -1;
            }

            // Lower triangle expand
            for
            (
                label indexC = losortCoarseStart[ir];
                indexC < losortCoarseStart[ir + 1];
                indexC++
            )
            {
                coeffLabel[lowerCoarseAddr[losortCoarseAddr[indexC]]] = -1;
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
    const labelList& rowR = crR.rowStart();
    const labelList& colR = crR.column();

    // Coefficients of restriction
    const scalarField& coeffR = R.coeffs();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    // Multiply the residual with restriction weights to obtain the initial
    // coarse residual

    for (label i = 0; i < nCoarseEqns_; i++)
    {
        for (label k = rowR[i]; k < rowR[i + 1]; k++)
        {
            // Multiply each coeff in row of restriction with the corresponding
            // residual coefficient (col index of R is the same as row index
            // of residual)

            // Col of R
            label j = colR[k];

            coarseRes[i] += mult(coeffR[k], res[j]);
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
    const labelList& rowP = crP.rowStart();
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
        for (label k = rowP[i]; k < rowP[i + 1]; k++)
        {
            // Multiply each coeff in row of prolongation with the corresponding
            // coarse residual coefficient (col index of prolongation must be
            // the same as row index of coarse residual coeff)

            // Col of P
            label j = colP[k];

            x[i] += mult(coeffP[k], coarseX[j]);
        }
    }
}


// ************************************************************************* //
