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

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    selectionAmgPolicy

Description
    Stueben-type selective AMG policy

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "selectionAmgPolicy.H"
#include "crMatrix.H"
#include "amgMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "SAMGInterfaceField.H"
#include "PriorityList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(selectionAmgPolicy, 0);

    addToRunTimeSelectionTable(amgPolicy, selectionAmgPolicy, matrix);
} // End namespace Foam


const Foam::debug::tolerancesSwitch Foam::selectionAmgPolicy::epsilon_
(
    "samgEpsilon",
    0.2
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::crMatrix> Foam::selectionAmgPolicy::filterProlongation
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


void Foam::selectionAmgPolicy::calcCoarsening()
{
//------------------------------------------------------------------------------
//             MATRIX DATA: ADDRESSING, COEFFICIENTS
//------------------------------------------------------------------------------

    // Get addressing
    const label nRows = matrix().lduAddr().size();
    const unallocLabelList& lowerAddr = matrix().lduAddr().lowerAddr();
    const unallocLabelList& upperAddr = matrix().lduAddr().upperAddr();
    const unallocLabelList& losortAddr = matrix().lduAddr().losortAddr();
    const unallocLabelList& ownerStart = matrix().lduAddr().ownerStartAddr();
    const unallocLabelList& losortStart = matrix().lduAddr().losortStartAddr();

    // Get matrix coefficients
    const scalarField& diag = matrix().diag();
    const scalarField& upper = matrix().upper();
    const scalarField& lower = matrix().lower();

//------------------------------------------------------------------------------
//        CRITERIA FOR COARSENING: FIND STRONG CONNECTIONS FOR EACH ROW
//------------------------------------------------------------------------------

    // STEP 1
    // Find the largest negative coefficient (sign opposite to diagonal
    // coefficient) in each row.
    // Multiply this largest coefficient with epsilon_ and compare all negative
    // coefficients in row to this value (epsilonStrongCoeff).
    // If they are larger than epsilonStrongCoeff, they are called
    // INFLUENCES of the equation, or STRONG connections.
    // If they are smaller than epsilonStrongCoeff, they are WEAK connections.
    // Store all strong connections into a matrix which will contain only the
    // strong negative connections for each row.

    // Field for storing the largest negative coefficient in each row
    scalarField epsilonStrongCoeff(nRows, 0);

    // Select the strongest coefficient in each row
    // (There is no need to operate in rows in the first 2 loops, but the
    // code is left row-wise for clarity)
    for (label i = 0; i < nRows; i++)
    {
        const scalar signDiag = sign(diag[i]);

        // Do lower triangle coefficient for the row
        for (label jp = losortStart[i]; jp < losortStart[i + 1]; jp++)
        {
            const scalar magAij = mag(min(signDiag*lower[losortAddr[jp]], 0));

            if (magAij > epsilonStrongCoeff[i])
            {
                epsilonStrongCoeff[i] = magAij;
            }
        }

        // Do upper triangle coefficient for the row
        for (label ip = ownerStart[i]; ip < ownerStart[i + 1]; ip++)
        {
            const scalar magAij = mag(min(signDiag*upper[ip], 0));

            if (magAij > epsilonStrongCoeff[i])
            {
                epsilonStrongCoeff[i] = magAij;
            }
        }
    }

    // Multiply the largest negative coefficient in row with epsilon_
    epsilonStrongCoeff *= epsilon_();

    // Count coefficients in row stronger than fraction of the strongest
    // (epsilonStrongCoeff_)
    labelList strongCoeffCounter(nRows, 0);

    for (label i = 0; i < nRows; i++)
    {
        const scalar signDiag = sign(diag[i]);

        // Do lower triangle coefficient for the row
        for (label jp = losortStart[i]; jp < losortStart[i + 1]; jp++)
        {
            const scalar magAij = mag(min(signDiag*lower[losortAddr[jp]], 0));

            if (magAij > epsilonStrongCoeff[i])
            {
                strongCoeffCounter[i]++;
            }
        }

        // Do upper triangle coefficient for the row
        for (label ip = ownerStart[i]; ip < ownerStart[i + 1]; ip++)
        {
            const scalar magAij = mag(min(signDiag*upper[ip], 0));

            if (magAij > epsilonStrongCoeff[i])
            {
                strongCoeffCounter[i]++;
            }
        }
    }

    // Create a crMatrix that will contain only the strong connections
    // (influences) of the equations: the matrix is called STRONG
    crMatrix strong(nRows, nRows, strongCoeffCounter);

    // Set addressing for the matrix:
    // stongCol and strongCoeff are arrays needed to create a compressed row
    // matrix
    const labelList& strongRow = strong.crAddr().rowStart();
    labelList& strongCol = strong.column();
    scalarField& strongCoeff = strong.coeffs();

    // Reset strongCoeffCounter for re-use
    strongCoeffCounter = 0;

    // Fill in the empty arrays of the strong matrix
    for (label i = 0; i < nRows; i++)
    {
        const scalar signDiag = sign(diag[i]);

        // Do lower triangle coefficient for the row
        for (label jp = losortStart[i]; jp < losortStart[i + 1]; jp++)
        {
            const scalar magAij = mag(min(signDiag*lower[losortAddr[jp]], 0));

            if (magAij > epsilonStrongCoeff[i])
            {
                strongCol[strongRow[i] + strongCoeffCounter[i]] =
                    lowerAddr[losortAddr[jp]];
                strongCoeff[strongRow[i] + strongCoeffCounter[i]] =
                    lower[losortAddr[jp]];

                strongCoeffCounter[i]++;
            }
        }

        // Do upper triangle coefficient for the row
        for (label ip = ownerStart[i]; ip < ownerStart[i + 1]; ip++)
        {
            const scalar magAij = mag(min(signDiag*upper[ip], 0));

            if (magAij > epsilonStrongCoeff[i])
            {
                strongCol[strongRow[i] + strongCoeffCounter[i]] = upperAddr[ip];
                strongCoeff[strongRow[i] + strongCoeffCounter[i]] = upper[ip];
                strongCoeffCounter[i]++;
            }
        }
    }

//------------------------------------------------------------------------------
//           COARSENING: SORT EQUATIONS INTO COARSE AND FINE SUBSETS
//------------------------------------------------------------------------------

    // STEP 2
    // Transpose the strong matrix: now columns of the strong matrix have become
    // rows of the transposed strong matrix.
    // The coefficients in the row of the transposed strong matrix are called
    // DEPENDENCIES of the equation.
    // A weight is assigned to each equation: the number of dependencies. All
    // equations are UNDECIDED at this point.
    // The equation with the largest weight is chosen to be coarse, all its
    // undecided neighbours in the transposed strong matrix are fine. The
    // neighbours of these fine equations have their weights incresed. Weights
    // of the neighbours of the coarse equation in the strong matrix are
    // decremented.
    // This results in the coarse-fine splitting.

    // Transpose the strong matrix to use for coarsening
    crAddressing Taddr = strong.crAddr().T();
    const labelList& tRow = Taddr.rowStart();
    const labelList& tCol = Taddr.column();

    // Create array to store the coarse/fine label. At this point, all equations
    // are undecided.
    rowLabel_.setSize(nRows, UNDECIDED);

    // Set initial equationWeight in the special type of list called
    // PriorityList.
    PriorityList<label> equationWeight(nRows);

    for (label i = 0; i < nRows; i++)
    {
        // Set weights for each equation: equations that my equation is strongly
        // influencing (dependencies).
        equationWeight.set(i, tRow[i + 1] - tRow[i]);
    }

    // Mark disconnected rows as fine. This removes solo equations which are now
    // solved only on the fine level.
    for (label i = 0; i < nRows; i++)
    {
        if (strongRow[i + 1] == strongRow[i])
        {
            rowLabel_[i] = FINE;
        }
    }

    // Count coarse equations
    nCoarseEqns_ = 0;

    while (!equationWeight.empty())
    {
        label topElement = equationWeight.removeHead();

        if (rowLabel_[topElement] == UNDECIDED)
        {
            // Make highest weight equation coarse
            rowLabel_[topElement] = nCoarseEqns_;
            nCoarseEqns_++;

            // Decrement weights of influences of my coarse equation. Since it
            // is coarse, it does not need interpolation from them, so a smaller
            // weight decreases their chance to become coarse.
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
                    // When weight is updated, the PriorityList
                    // will re-sort its contents
                    equationWeight.updateWeight
                    (
                        j,
                        equationWeight.weights()[j] - 1
                    );
                }
            }

            // All undecided dependencies of the coarse equation become fine:
            // they will interpolate from the coarse equation.
            // The interior for loop goes through their influences: if they are
            // undecided, their weight is incremented. They have more chance to
            // become coarse, thus our fine equation would have more coarse
            // influences to interpolate from.
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
                            // Note: When weight is updated, the PriorityList
                            // will re-sort its contents
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

    // STEP 3
    // The equations are split into coarse and fine which is stored in rowLabel_
    // list: fine equations are marked with -1, coarse equations are marked with
    // the index of the coarse equation (0, 1, 2,...).
    // For calculating the interpolation weights, it is neccesary to take into
    // account strong fine neighbours: it is done by summing all strong
    // connections (num) and dividing them by the sum of strong coarse
    // connections (den). This is the scaling factor (alpha) which is always
    // larger or equal to 1.
    // In this part, we will only take the coefficients on the coupled processor
    // boundary. Internal coefficents from strong matrix are considered in the
    // next step.

    // Sum of strong fine and coarse connections
    scalarField num(nRows, 0);

    // Sum of weak connections (positive coefficients) in row (sign equal to
    // diagonal)
    scalarField Dii(sign(diag)*diag);

    scalar Dij, den, signDii;

    // Coupled boundary treatment
    // To take into account the neighbours across processor boundary, include
    // them into the scaling factor.
    // Important:
    // bouCoeffs() = upper triangular coupled coefficient
    // intCoeffs() = lower triangular coupled coefficient
    // The sign of bouCoeffs and intCoeffs is opposite of the sign of internal
    // off-diag coeffs (they are on rhs).
    // On the local processor, we only account for upper coefficients, as the
    // column coeff is handled by the other processor.

    forAll (interfaceFields(), intI)
    {
        if (interfaceFields().set(intI))
        {
            // Get bouCoeffs
            const scalarField& cplUpper = bouCoeffs()[intI];

            // Get addressing
            const labelList& faceCells =
                interfaceFields()[intI].coupledInterface().faceCells();

            forAll (cplUpper, coeffI)
            {
                // Get row index from faceCells
                const label i = faceCells[coeffI];

                // Handle the unknown sign of diagonal: multiplying
                // row coeffs by the sign
                signDii = sign(diag[i]);

                // Adjust sign of off-diag coeff
                // Note: additional minus because the sign of interface
                // coeffs is opposite from the normal matrix off-diagonal
                Dij = -signDii*cplUpper[coeffI];

                // Add negative coeff to num
                // Add positive coeff to diag to eliminate it
                num[i] += Foam::min(Dij, 0);
                Dii[i] += Foam::max(Dij, 0);
            }
        }
    }

//------------------------------------------------------------------------------
//                    CALCULATING PROLONGATION WEIGHTS
//------------------------------------------------------------------------------

    // STEP 4
    // Calculate the interpolation weight, i.e. how much my fine equation gets
    // from each of its strong coarse neighbours.
    // The formula for interpolation:
    // w = - alpha * Dii^-1 * Aij_strongCoarse
    // In this part, calculation of alpha (num/den) is finished by looping
    // through the fine matrix (num) and strong matrix (den).

    // Calculate the maximal number of interpolation weights that can appear:
    // one for each strong neighbour of fine equations, and only one for coarse
    // equations (injection from a single equation).
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

    // Lists for assembling the prolongation matrix
    labelList pRow(nRows + 1);
    labelList pCol(maxPCount, 0);
    scalarField pCoeff(maxPCount, 0);

    // Start row assembly
    pRow[0] = 0;

    for (label i = 0; i < nRows; i++)
    {
        label rowCount = pRow[i];

        // Calculate the scaling factor alpha (num/den)

        // Handle the unknown sign of diagonal: multiplying row coeffs
        // by the sign of the diagonal.
        signDii = sign(diag[i]);

        // Do lower triangle coefficient for the row
        for (label jp = losortStart[i]; jp < losortStart[i + 1]; jp++)
        {
            // Adjust sign of off-diag coeff
            Dij = signDii*lower[losortAddr[jp]];

            // Add negative coeff to num
            // Add positive coeff to diag to eliminate it
            num[i] += Foam::min(Dij, 0);
            Dii[i] += Foam::max(Dij, 0);
        }

        // Do upper triangle coefficient for the row
        for (label ip = ownerStart[i]; ip < ownerStart[i + 1]; ip++)
        {
            // Adjust sign of off-diag coeff
            Dij = signDii*upper[ip];

            // Add negative coeff to num
            // Add positive coeff to diag to eliminate it
            num[i] += Foam::min(Dij, 0);
            Dii[i] += Foam::max(Dij, 0);
        }

        // Completed calculating numerator for the scaling factor.
        // Calculate denominator and weights.
        if (rowLabel_[i] == FINE)
        {
            // Calculate denominator of scaling factor.
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

            // Calculate interpolation weights.
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
                    // Prolongation coefficient: this formula preserves the fine
                    // matrix row sums. This means that for:
                    // - diagonally equal equations, the sum of weights in the
                    // corresponding row of the prolongation matrix will be 1.
                    // - diagonally dominant equations, the sum of weights in
                    // the corresponding row of prolongation matrix will be < 1.
                    // - diagonally weak equations, the sum of weights in the
                    // corresponding row of prolongation matrix will be > 1.
                    pCoeff[rowCount] = -(num[i]/den)*strongCoeff[sip]/Dii[i];

                    // Prolongation coefficent scaled for rows to sum into 1:
                    // we have to assume all the equations are diagonally equal
                    // in order to get the row of prolongation to sum into 1.
                    // num[i] is a sum of all negative (sign opposite to
                    // diagonal) off-diagonal coefficients and SHOULD be equal
                    // to the diagonal. That is why -num[i] and Dii[i] have
                    // disappeared. Note: if there are positive connections in
                    // the row, this is not valid.

                    //pCoeff[rowCount] = strongCoeff[sip]/den;

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

    // Copy coefficients
    prolongation.coeffs().transfer(pCoeff);

    // Check prolongation matrix
    if (lduMatrix::debug > 2)
    {
        label nSoloEqn = 0;

        scalarField sumRow(nRows, 0);
        const labelList& prolongationRow = Pptr_->crAddr().rowStart();
        const scalarField& prolongationCoeff = Pptr_->coeffs();

        for (label rowI = 0; rowI < nRows; rowI++)
        {
            if (prolongationRow[rowI] == prolongationRow[rowI + 1])
            {
                Pout<< "Solo eqn: " << rowI << " original "
                    << ownerStart[rowI + 1] - ownerStart[rowI] << " "
                    << losortStart[rowI + 1] - losortStart[rowI] << " C/F: "
                    << rowLabel_[rowI] << " strong "
                    << strongRow[rowI + 1] - strongRow[rowI]
                    << endl;

                nSoloEqn++;
            }

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

        Pout<< "sumRow (min, max) = (" << min(sumRow) << " "
            << max(sumRow) << ").  nSolo = " << nSoloEqn
            << " weight (min, max) = (" << min(prolongationCoeff)
            << " " << max(prolongationCoeff) << ")"
            << endl;
    }


    // The decision on parallel agglomeration needs to be made for the
    // whole gang of processes; otherwise I may end up with a different
    // number of agglomeration levels on different processors.
    if
    (
        nCoarseEqns_ > minCoarseEqns()
     && 3*nCoarseEqns_ <= 2*nRows
    )
    {
        coarsen_ = true;
    }

    reduce(coarsen_, andOp<bool>());

    if (lduMatrix::debug >= 2)
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
        Rptr_ = new crMatrix(Pptr_->T());
    }
    else
    {
        // Coarsening did not succeed.  Delete Pptr
        deleteDemandDrivenData(Pptr_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::selectionAmgPolicy::selectionAmgPolicy
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& bouCoeffs,
    const FieldField<Field, scalar>& intCoeffs,
    const lduInterfaceFieldPtrsList& interfaceFields,
    const label groupSize,
    const label minCoarseEqns
)
:
    amgPolicy
    (
        matrix,
        bouCoeffs,
        intCoeffs,
        interfaceFields,
        groupSize,
        minCoarseEqns
    ),
    nCoarseEqns_(0),
    coarsen_(false),
    Pptr_(nullptr),
    Rptr_(nullptr)
{
    calcCoarsening();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::selectionAmgPolicy::~selectionAmgPolicy()
{
    deleteDemandDrivenData(Rptr_);
    deleteDemandDrivenData(Pptr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::amgMatrix> Foam::selectionAmgPolicy::restrictMatrix() const
{
    if (!coarsen_)
    {
        FatalErrorIn
        (
            "autoPtr<amgMatrix> selectionAmgPolicy::restrictMatrix() const"
        )   << "Requesting coarse matrix when it cannot be created"
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
    const label nEqns = matrix().lduAddr().size();

    if
    (
        crR.nCols() != nEqns
     || crP.nRows() != nEqns
     || crR.nRows() != nCoarseEqns_
    )
    {
        FatalErrorIn
        (
            "autoPtr<amgMatrix> Foam::selectionAmgPolicy::restrictMatrix"
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
    const unallocLabelList& rowA = matrix().lduAddr().ownerStartAddr();
    const unallocLabelList& upperAddr = matrix().lduAddr().upperAddr();

    // Addressing for lower triangle loop
    const unallocLabelList& lowerAddr = matrix().lduAddr().lowerAddr();
    const unallocLabelList& losortAddr = matrix().lduAddr().losortAddr();
    const unallocLabelList& losortStart = matrix().lduAddr().losortStartAddr();

    // Prolongation addressing
    const labelList& rowP = crP.rowStart();
    const labelList& colP = crP.column();

    // In order to avoid searching for the off-diagonal coefficient,
    // a mark array is used for each row's assembly.
    // Mark records the index of the off-diagonal coefficient in the row
    // for each neighbour entry. It is reset after completing each row of the
    // coarse matrix.
    labelList coeffLabel(nCoarseEqns_, -1);

    // Create coarse addressing: record neighbours for every row
    // Note: neighbour indices are added as they appear and will be sorted
    // on completion.
    List<labelHashSet> coarseNbrsSets(nCoarseEqns_);

//------------------------------------------------------------------------------
//     COUNT COARSE COEFFICIENTS IN TRIPLE PRODUCT AND CREATE ADDRESSING
//------------------------------------------------------------------------------

    // This part of the code is specialised for symmetric and asymmetric matrix
    // The code is identical, but, for symmetric matrix, we do not store the
    // contributions to the lower triangle

    // Symmetric matrix
    if (matrix().symmetric())
    {
        // This loop will be used to count the number of coeffs appearing in the
        // COARSE matrix and creating the owner-neighbour addressing of the
        // coarse matrix.

        // Dynamic Lists are used to create coarseOwner and coarseNeighbour
        // arrays in which owner and neighbour are stored for the upper tringle
        // Note: the neighbour array should be sorted before constructing the
        // coarseMatrix

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
                            // SKIP!
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
                            // SKIP!
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
                        // SKIP!
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
    }
    // Asymmetric matrix
    else
    {
        // This loop will be used to count the number of coeffs appearing in the
        // COARSE matrix and creating the owner-neighbour addressing of the
        // coarse matrix.

        // Dynamic Lists are used to create coarseOwner and coarseNeighbour
        // arrays in which owner and neighbour are stored for the upper tringle
        // Note: the neighbour array should be sorted before constructing the
        // coarseMatrix

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
    lduPrimitiveMesh* coarseAddrPtr =
        new lduPrimitiveMesh
        (
            nCoarseEqns_,
            coarseOwner,
            coarseNeighbour,
            false // Do not reuse storage
        );


//------------------------------------------------------------------------------
//                       CREATE COARSE MATRIX INTERFACES
//------------------------------------------------------------------------------

    // Create coarse-level coupled interface fields

    // Set the coarse interfaces and coefficients
    lduInterfacePtrsList* coarseInterfacesPtr =
        new lduInterfacePtrsList(interfaceFields().size());
    lduInterfacePtrsList& coarseInterfaces = *coarseInterfacesPtr;

    // Set the coarse interfaceFields and coefficients
    lduInterfaceFieldPtrsList* coarseInterfaceFieldsPtr =
        new lduInterfaceFieldPtrsList(interfaceFields().size());
    lduInterfaceFieldPtrsList& coarseInterfaceFields =
        *coarseInterfaceFieldsPtr;

    FieldField<Field, scalar>* coarseBouCoeffsPtr =
        new FieldField<Field, scalar>(interfaceFields().size());
    FieldField<Field, scalar>& coarseBouCoeffs = *coarseBouCoeffsPtr;

    FieldField<Field, scalar>* coarseIntCoeffsPtr =
        new FieldField<Field, scalar>(interfaceFields().size());
    FieldField<Field, scalar>& coarseIntCoeffs = *coarseIntCoeffsPtr;

    labelListList coarseInterfaceAddr(interfaceFields().size());

    // Store owner and neighbour prolongation to avoid tangled communications
    PtrList<crMatrix> ownInterfaceProlongation(interfaceFields().size());
    PtrList<crMatrix> nbrInterfaceProlongation(interfaceFields().size());

    // Initialise transfer of matrix prolongation on the interface
    forAll (interfaceFields(), intI)
    {
        if (interfaceFields().set(intI))
        {
            const labelList& fineFaceCells =
                interfaceFields()[intI].coupledInterface().faceCells();

            // Filter local prolongation matrix and return to
            // filteredProlongation
            ownInterfaceProlongation.set
            (
                intI,
                filterProlongation(P, fineFaceCells).ptr()
            );

            // Send the filtered prolongation matrix across
            interfaceFields()[intI].coupledInterface().initProlongationTransfer
            (
                Pstream::blocking,
                ownInterfaceProlongation[intI]
            );
        }
    }

    forAll (interfaceFields(), intI)
    {
        if (interfaceFields().set(intI))
        {
            const lduInterface& fineInterface =
                interfaceFields()[intI].coupledInterface();

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

    // Add the coarse level SAMG interfaces
    forAll (interfaceFields(), intI)
    {
        if (interfaceFields().set(intI))
        {
            const lduInterface& fineInterface =
                interfaceFields()[intI].coupledInterface();

            // Note: select SAMG interface for matrix selection
            // Use filtered prolongation on master and slave side
            coarseInterfaces.set
            (
                intI,
                SAMGInterface::New
                (
                    *coarseAddrPtr,
                    ownInterfaceProlongation[intI],
                    coarseInterfaces,
                    fineInterface,
                    nbrInterfaceProlongation[intI]
                ).ptr()
            );
        }
    }

    forAll (interfaceFields(), intI)
    {
        if (interfaceFields().set(intI))
        {
            const SAMGInterface& coarseInterface =
                refCast<const SAMGInterface>(coarseInterfaces[intI]);

            // Note: select SAMG interface field for matrix selection
            coarseInterfaceFields.set
            (
                intI,
                SAMGInterfaceField::New
                (
                    coarseInterface,
                    interfaceFields()[intI]
                ).ptr()
            );

            // Note: scalar agglomeration is done by the interface
            // (always scalar) but in the block matrix it is done by a
            // templated block interface field
            coarseBouCoeffs.set
            (
                intI,
                coarseInterface.selectCoeffs(bouCoeffs()[intI])
            );

            coarseIntCoeffs.set
            (
                intI,
                coarseInterface.selectCoeffs(intCoeffs()[intI])
            );

            coarseInterfaceAddr[intI] = coarseInterface.faceCells();
        }
    }

    // Add interfaces
    coarseAddrPtr->addInterfaces
    (
        *coarseInterfacesPtr,
        coarseInterfaceAddr,
        matrix().patchSchedule()
    );

//------------------------------------------------------------------------------
//                            CREATE COARSE MATRIX
//------------------------------------------------------------------------------

    // Set the coarse level matrix
    lduMatrix* coarseMatrixPtr = new lduMatrix(*coarseAddrPtr);
    lduMatrix& coarseMatrix = *coarseMatrixPtr;

    // This part of the code is specialised for symmetric and asymmetric matrix
    // The code is identical, but, for symmetric matrix, we do not store the
    // contributions to the lower triangle

    // Symmetric matrix
    if (matrix().symmetric())
    {
        scalarField& coarseUpper = coarseMatrix.upper();
        scalarField& coarseDiag = coarseMatrix.diag();

//------------------------------------------------------------------------------
//                           GET COEFFICIENTS
//------------------------------------------------------------------------------

        // Coefficients of matrix A
        const scalarField& diag = matrix().diag();
        const scalarField& upper = matrix().upper();
        const scalarField& lower = matrix().lower();

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

        // Get coarse matrix addressing
        const unallocLabelList& rowC = coarseMatrix.lduAddr().ownerStartAddr();
        const unallocLabelList& upperCoarseAddr =
            coarseMatrix.lduAddr().upperAddr();

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
                    const scalar ra = lower[faceA]*coeffR[indexR];

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
                        }
                        else if (ir == jp)
                        {
                            // Found COARSE diagonal coefficient
                            coarseDiag[ir] += ra*coeffP[indexP];
                        }
                        else
                        {
                            // Found upper COARSE triangle
                            label face = coeffLabel[jp];
                            coarseUpper[face] += ra*coeffP[indexP];
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

                    const scalar ra = coeffR[indexR]*upper[faceA];

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
                        }
                        else if (ir == jp)
                        {
                           // Found COARSE diagonal
                            coarseDiag[ir] += ra*coeffP[indexP];
                        }
                        else
                        {
                            // Found upper COARSE triangle
                            label face = coeffLabel[jp];
                            coarseUpper[face] += ra*coeffP[indexP];
                        }
                    }
                }

                // DIAGONAL
                const scalar ra = coeffR[indexR]*diag[jr];

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
                    }
                    else if (ir == jp)
                    {
                        // Found COARSE diagonal
                        coarseDiag[ir] += ra*coeffP[indexP];
                    }
                    else
                    {
                        // Found upper COARSE triangle
                        label face = coeffLabel[jp];
                        coarseUpper[face] += ra*coeffP[indexP];
                    }
                }
            }

            // Clean out coeffLabel
            // Upper triangle expand
            for (label coarseK = rowC[ir]; coarseK < rowC[ir + 1]; coarseK++)
            {
                coeffLabel[upperCoarseAddr[coarseK]] = -1;
            }
        }
    }
    // Asymmetric matrix
    else
    {
        scalarField& coarseUpper = coarseMatrix.upper();
        scalarField& coarseDiag = coarseMatrix.diag();
        scalarField& coarseLower = coarseMatrix.lower();

//------------------------------------------------------------------------------
//                           GET COEFFICIENTS
//------------------------------------------------------------------------------

        // Coefficients of matrix A
        const scalarField& diag = matrix().diag();
        const scalarField& upper = matrix().upper();
        const scalarField& lower = matrix().lower();

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
                    const scalar ra = lower[faceA]*coeffR[indexR];

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
                            coarseLower[face] += ra*coeffP[indexP];
                        }
                        else if (ir == jp)
                        {
                            // Found COARSE diagonal coefficient
                            coarseDiag[ir] += ra*coeffP[indexP];
                        }
                        else
                        {
                            // Found upper COARSE triangle
                            label face = coeffLabel[jp];
                            coarseUpper[face] += ra*coeffP[indexP];
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

                    const scalar ra = coeffR[indexR]*upper[faceA];

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
                            coarseLower[face] += ra*coeffP[indexP];
                        }
                        else if (ir == jp)
                        {
                           // Found COARSE diagonal
                            coarseDiag[ir] += ra*coeffP[indexP];
                        }
                        else
                        {
                            // Found upper COARSE triangle
                            label face = coeffLabel[jp];
                            coarseUpper[face] += ra*coeffP[indexP];
                        }
                    }
                }

                // DIAGONAL
                const scalar ra = coeffR[indexR]*diag[jr];

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
                        coarseLower[face] += ra*coeffP[indexP];
                    }
                    else if (ir == jp)
                    {
                        // Found COARSE diagonal
                        coarseDiag[ir] += ra*coeffP[indexP];
                    }
                    else
                    {
                        // Found upper COARSE triangle
                        label face = coeffLabel[jp];
                        coarseUpper[face] += ra*coeffP[indexP];
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
    }

    // Create and return coarse level
    return autoPtr<amgMatrix>
    (
        new amgMatrix
        (
            coarseAddrPtr,
            coarseInterfacesPtr,
            coarseMatrixPtr,
            coarseBouCoeffsPtr,
            coarseIntCoeffsPtr,
            coarseInterfaceFieldsPtr
        )
    );
}


void Foam::selectionAmgPolicy::restrictResidual
(
    const scalarField& res,
    scalarField& coarseRes
) const
{
    coarseRes = 0;
    Rptr_->dotPlus(coarseRes, res);
}


void Foam::selectionAmgPolicy::prolongateCorrection
(
    scalarField& x,
    const scalarField& coarseX
) const
{
    // No need to initialise x into 0 because we are correcting the solution
    // based on the error calculated from residual eqn Ae = r.
    // Correct previously obtained solution with the correction term from
    // prolongating the coarse level error.
    Pptr_->dotPlus(x, coarseX);
}


// ************************************************************************* //
