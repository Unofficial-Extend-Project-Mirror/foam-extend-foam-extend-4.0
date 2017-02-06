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
    BlockMatrixAgglomeration

Description
    Agglomerative block matrix AMG corsening

Author
    Klas Jareteg, 2012-12-13

\*---------------------------------------------------------------------------*/

#include "BlockMatrixAgglomeration.H"
#include "boolList.H"
#include "tolerancesSwitch.H"
#include "coeffFields.H"
#include "addToRunTimeSelectionTable.H"
#include "BlockAMGInterfaceField.H"
#include "coarseBlockAMGLevel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::debug::tolerancesSwitch
Foam::BlockMatrixAgglomeration<Type>::weightFactor_
(
    "aamgWeightFactor",
    0.65
);


template<class Type>
const Foam::debug::tolerancesSwitch
Foam::BlockMatrixAgglomeration<Type>::diagFactor_
(
    "aamgDiagFactor",
    1e-8
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::BlockMatrixAgglomeration<Type>::calcAgglomeration()
{
    // Algorithm:
    // 1) Create temporary equation addressing using a double-pass algorithm.
    //    to create the offset table.
    // 2) Loop through all equations and for each equation find the best fit
    //    neighbour.  If all neighbours are grouped, add equation to best group

    // Create row-based addressing

    const label nRows = matrix_.lduAddr().size();

    const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();

    // For each equation get number of coefficients in a row
    labelList cols(upperAddr.size() + lowerAddr.size());
    labelList cIndex(upperAddr.size() + lowerAddr.size());
    labelList rowOffset(nRows + 1, 0);

    // Count the number of coefficients
    forAll (upperAddr, coeffI)
    {
        rowOffset[upperAddr[coeffI]]++;
    }

    forAll (lowerAddr, coeffI)
    {
        rowOffset[lowerAddr[coeffI]]++;
    }

    label nCoeffs = 0;

    forAll (rowOffset, eqnI)
    {
        nCoeffs += rowOffset[eqnI];
    }

    rowOffset[nRows] = nCoeffs;

    for (label eqnI = nRows - 1; eqnI >= 0; --eqnI)
    {
        rowOffset[eqnI] = rowOffset[eqnI  + 1] - rowOffset[eqnI];
    }

    rowOffset[0] = 0;

    // Create column and coefficient index array
    {
        // Use agglomIndex to count number of entries per row.
        // Reset the list for counting
        labelList& nPerRow = agglomIndex_;
        nPerRow = 0;

        forAll (upperAddr, coeffI)
        {
            cols[rowOffset[upperAddr[coeffI]] + nPerRow[upperAddr[coeffI]]] =
                lowerAddr[coeffI];

            cIndex[rowOffset[upperAddr[coeffI]] + nPerRow[upperAddr[coeffI]]] =
                coeffI;

            nPerRow[upperAddr[coeffI]]++;
        }

        forAll (lowerAddr, coeffI)
        {
            cols[rowOffset[lowerAddr[coeffI]] + nPerRow[lowerAddr[coeffI]]] =
                upperAddr[coeffI];

            cIndex[rowOffset[lowerAddr[coeffI]] + nPerRow[lowerAddr[coeffI]]] =
                coeffI;

            nPerRow[lowerAddr[coeffI]]++;
        }

        // Reset agglomeration index array
        agglomIndex_ = -1;
    }


    // Calculate agglomeration

    // Get matrix coefficients
    const CoeffField<Type>& diag = matrix_.diag();

    // Coefficient magnitudes are pre-calculated
    scalarField magDiag(diag.size());
    scalarField magOffDiag(upperAddr.size());

    normPtr_->coeffMag(diag, magDiag);

    if (matrix_.asymmetric())
    {
        scalarField magUpper(upperAddr.size());
        scalarField magLower(upperAddr.size());

        normPtr_->coeffMag(matrix_.upper(), magUpper);
        normPtr_->coeffMag(matrix_.lower(), magLower);

        magOffDiag = Foam::max(magUpper, magLower);
    }
    else if (matrix_.symmetric())
    {
        normPtr_->coeffMag(matrix_.upper(), magOffDiag);
    }
    else
    {
        // Diag only matrix.  Reset and return
        agglomIndex_ = 0;
        nCoarseEqns_ = 1;

        return;
    }

    labelList sizeOfGroups(nRows, 0);

    nCoarseEqns_ = 0;

    // Gather disconnected and weakly connected equations into cluster zero
    // Weak connection is assumed to be the one where the off-diagonal
    // coefficient is smaller than diagFactor_*diag
    {
        // Algorithm
        // Mark all cells to belong to zero cluster
        // Go through all upper and lower coefficients and for the ones
        // larger than threshold mark the equations out of cluster zero

        scalarField magScaledDiag = diagFactor_()*magDiag;

        boolList zeroCluster(diag.size(), true);

        forAll (magOffDiag, coeffI)
        {
            if (magOffDiag[coeffI] > magScaledDiag[upperAddr[coeffI]])
            {
                zeroCluster[upperAddr[coeffI]] = false;
            }

            if (magOffDiag[coeffI] > magScaledDiag[lowerAddr[coeffI]])
            {
                zeroCluster[lowerAddr[coeffI]] = false;
            }
        }

        // Collect solo equations
        forAll (zeroCluster, eqnI)
        {
            if (zeroCluster[eqnI])
            {
                // Found solo equation
                nSolo_++;

                agglomIndex_[eqnI] = nCoarseEqns_;
            }
        }

        if (nSolo_ > 0)
        {
            // Found solo equations
            nCoarseEqns_++;

            if (BlockLduMatrix<Type>::debug >= 3)
            {
                Pout<< "Found " << nSolo_ << " weakly connected equations."
                    << endl;
            }
        }
    }

    // Go through the off-diagonal and create clusters, marking the child array
    label indexUngrouped, indexGrouped;
    label colI, curEqn, nextEqn, groupPassI;

    scalar magRowDiag, magColDiag;
    scalar weight, weightUngrouped, weightGrouped;

    for (label eqnI = 0; eqnI < nRows; eqnI++)
    {
        if (agglomIndex_[eqnI] == -1)
        {
            curEqn = eqnI;

            indexUngrouped = -1;
            indexGrouped = -1;

            agglomIndex_[curEqn] = nCoarseEqns_;

            magRowDiag = magDiag[curEqn];

            for (groupPassI = 1; groupPassI < groupSize_; groupPassI++)
            {
                weightUngrouped = 0;
                weightGrouped = 0;

                indexUngrouped = -1;
                indexGrouped = -1;

                for
                (
                    label rowCoeffI = rowOffset[curEqn];
                    rowCoeffI < rowOffset[curEqn + 1];
                    rowCoeffI++
                )
                {
                    colI = cols[rowCoeffI];

                    magColDiag = magDiag[colI];

                    weight = magOffDiag[cIndex[rowCoeffI]]/
                        max(magRowDiag, magColDiag);

                    if (agglomIndex_[colI] == -1)
                    {
                        if (indexUngrouped == -1 || weight > weightUngrouped)
                        {
                            indexUngrouped = rowCoeffI;
                            weightUngrouped = weight;
                        }
                    }
                    else if (agglomIndex_[curEqn] != agglomIndex_[colI])
                    {
                        if (indexGrouped == -1 || weight > weightGrouped)
                        {
                            indexGrouped = rowCoeffI;
                            weightGrouped = weight;
                        }
                    }
                }

                if
                (
                    indexUngrouped != -1
                 && (
                        indexGrouped == -1
                     || weightUngrouped >= weightFactor_()*weightGrouped
                    )
                )
                {
                    // Found new element of group.  Add it and use as
                    // start of next search

                    nextEqn = cols[indexUngrouped];

                    agglomIndex_[nextEqn] = agglomIndex_[curEqn];
                    sizeOfGroups[agglomIndex_[curEqn]]++;

                    curEqn = nextEqn;
                }
                else
                {
                    // Group full or cannot be extended
                    break;
                }
            }

            if
            (
                groupPassI > 1
             || indexGrouped == -1
             || (
                    sizeOfGroups[agglomIndex_[cols[indexGrouped]]]
                  > (groupSize_ + 2)
                )
            )
            {
                sizeOfGroups[agglomIndex_[eqnI]]++;
                nCoarseEqns_++;
            }
            else
            {
                agglomIndex_[eqnI] = agglomIndex_[cols[indexGrouped]];
                sizeOfGroups[agglomIndex_[cols[indexGrouped]]]++;
            }
        }
    }

    // The decision on parallel agglomeration needs to be made for the
    // whole gang of processes; otherwise I may end up with a different
    // number of agglomeration levels on different processors.

    // If the number of coarse equations is les than minimum and
    // if the matrix has reduced in size by at least 1/3, coarsen
    if
    (
        nCoarseEqns_ > BlockMatrixCoarsening<Type>::minCoarseEqns()
     && 3*nCoarseEqns_ <= 2*nRows
    )
    {
        coarsen_ = true;
    }

    reduce(coarsen_, andOp<bool>());

    if (BlockLduMatrix<Type>::debug >= 3)
    {
        Pout << "Coarse level size: " << nCoarseEqns_;

        if (coarsen_)
        {
            Pout << ".  Accepted" << endl;
        }
        else
        {
            Pout << ".  Rejected" << endl;
        }

        // Count cluster size
        labelList clusterSize(nCoarseEqns_, 0);

        forAll (agglomIndex_, eqnI)
        {
            clusterSize[agglomIndex_[eqnI]]++;
        }

        label minClusterSize = gMin(clusterSize);
        label maxClusterSize = gMax(clusterSize);

        Info<< "Cluster size: min = " << minClusterSize
            << " max = " << maxClusterSize << endl;
    }
}


template<class Type>
void Foam::BlockMatrixAgglomeration<Type>::restrictDiag
(
    const CoeffField<Type>& Coeff,
    CoeffField<Type>& coarseCoeff
) const
{
    typedef CoeffField<Type> TypeCoeffField;

    if
    (
        Coeff.activeType() == blockCoeffBase::SQUARE
     && coarseCoeff.activeType() == blockCoeffBase::SQUARE
    )
    {
        typedef typename TypeCoeffField::squareType squareType;
        typedef typename TypeCoeffField::squareTypeField squareTypeField;

        squareTypeField& activeCoarseCoeff = coarseCoeff.asSquare();
        const squareTypeField& activeCoeff = Coeff.asSquare();

        forAll (coarseCoeff, i)
        {
            activeCoarseCoeff[i] = pTraits<squareType>::zero;
        }

        forAll (Coeff, i)
        {
            activeCoarseCoeff[agglomIndex_[i]] += activeCoeff[i];
        }
    }
    else if
    (
        Coeff.activeType() == blockCoeffBase::LINEAR
     && coarseCoeff.activeType() == blockCoeffBase::LINEAR
    )
    {
        typedef typename TypeCoeffField::linearType linearType;
        typedef typename TypeCoeffField::linearTypeField linearTypeField;

        linearTypeField& activeCoarseCoeff = coarseCoeff.asLinear();
        const linearTypeField& activeCoeff = Coeff.asLinear();

        forAll (coarseCoeff, i)
        {
            activeCoarseCoeff[i] = pTraits<linearType>::zero;
        }

        forAll (Coeff, i)
        {
            activeCoarseCoeff[agglomIndex_[i]] += activeCoeff[i];
        }
    }
    else if
    (
        Coeff.activeType() == blockCoeffBase::SCALAR
     && coarseCoeff.activeType() == blockCoeffBase::SCALAR
    )
    {
        typedef typename TypeCoeffField::scalarType scalarType;
        typedef typename TypeCoeffField::scalarTypeField scalarTypeField;

        scalarTypeField& activeCoarseCoeff = coarseCoeff.asScalar();
        const scalarTypeField& activeCoeff = Coeff.asScalar();

        forAll (coarseCoeff, i)
        {
            activeCoarseCoeff[i] = pTraits<scalarType>::zero;
        }

        forAll (Coeff, i)
        {
            activeCoarseCoeff[agglomIndex_[i]] += activeCoeff[i];
        }
    }
    else
    {
        FatalErrorIn
        (
            "void  BlockMatrixAgglomeration<Type>::restrictDiag() const"
        )   << "Problem in coeff type morphing"
            << abort(FatalError);
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockMatrixAgglomeration<Type>::agglomerateCoeffs
(
    const labelList& coeffRestrictAddr,
    Field<DiagType>& activeCoarseDiag,
    Field<ULType>& activeCoarseUpper,
    const Field<ULType>& activeFineUpper,
    const Field<ULType>& activeFineUpperTranspose
) const
{
    // Does the matrix have solo equations
    bool soloEqns = nSolo_ > 0;

    // Get addressing
    const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();

    forAll(coeffRestrictAddr, fineCoeffI)
    {
        label rmUpperAddr = agglomIndex_[upperAddr[fineCoeffI]];
        label rmLowerAddr = agglomIndex_[lowerAddr[fineCoeffI]];

        // If the coefficient touches block zero and
        //  solo equations are present, skip it
        if (soloEqns && (rmUpperAddr == 0 || rmLowerAddr == 0))
        {
            continue;
        }

        label cCoeff = coeffRestrictAddr[fineCoeffI];

        if (cCoeff >= 0)
        {
            activeCoarseUpper[cCoeff] += activeFineUpper[fineCoeffI];
        }
        else
        {
            // Add the fine face coefficient into the diagonal
            // Note: upper and lower coeffs are transpose of
            // each other.  HJ, 28/May/2014
            activeCoarseDiag[-1 - cCoeff] +=
                activeFineUpper[fineCoeffI]
              + activeFineUpperTranspose[fineCoeffI];
        }
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockMatrixAgglomeration<Type>::agglomerateCoeffs
(
    const labelList& coeffRestrictAddr,
    Field<DiagType>& activeCoarseDiag,
    Field<ULType>& activeCoarseUpper,
    const Field<ULType>& activeFineUpper,
    Field<ULType>& activeCoarseLower,
    const Field<ULType>& activeFineLower
) const
{
    // Does the matrix have solo equations
    bool soloEqns = nSolo_ > 0;

    // Get addressing
    const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();

    forAll(coeffRestrictAddr, fineCoeffI)
    {
        label rmUpperAddr = agglomIndex_[upperAddr[fineCoeffI]];
        label rmLowerAddr = agglomIndex_[lowerAddr[fineCoeffI]];

        // If the coefficient touches block zero and
        //  solo equations are present, skip it
        if (soloEqns && (rmUpperAddr == 0 || rmLowerAddr == 0))
        {
            continue;
        }

        label cCoeff = coeffRestrictAddr[fineCoeffI];

        if (cCoeff >= 0)
        {
            activeCoarseUpper[cCoeff] += activeFineUpper[fineCoeffI];
            activeCoarseLower[cCoeff] += activeFineLower[fineCoeffI];
        }
        else
        {
            // Add the fine face coefficients into the diagonal.
            activeCoarseDiag[-1 - cCoeff] +=
                activeFineUpper[fineCoeffI]
              + activeFineLower[fineCoeffI];
        }
    }
}


template<class Type>
void Foam::BlockMatrixAgglomeration<Type>::restrictDiagDecoupled
(
    const CoeffField<Type>& Coeff,
    CoeffField<Type>& coarseCoeff
) const
{
    typedef CoeffField<Type> TypeCoeffField;

    if
    (
        Coeff.activeType() == blockCoeffBase::LINEAR
     && coarseCoeff.activeType() == blockCoeffBase::LINEAR
    )
    {
        typedef typename TypeCoeffField::linearType linearType;
        typedef typename TypeCoeffField::linearTypeField linearTypeField;

        linearTypeField& activeCoarseCoeff = coarseCoeff.asLinear();
        const linearTypeField& activeCoeff = Coeff.asLinear();

        forAll (coarseCoeff, i)
        {
            activeCoarseCoeff[i] = pTraits<linearType>::zero;
        }

        forAll (Coeff, i)
        {
            activeCoarseCoeff[agglomIndex_[i]] += activeCoeff[i];
        }
    }
    else if
    (
        Coeff.activeType() == blockCoeffBase::SCALAR
     && coarseCoeff.activeType() == blockCoeffBase::SCALAR
    )
    {
        typedef typename TypeCoeffField::scalarType scalarType;
        typedef typename TypeCoeffField::scalarTypeField scalarTypeField;

        scalarTypeField& activeCoarseCoeff = coarseCoeff.asScalar();
        const scalarTypeField& activeCoeff = Coeff.asScalar();

        forAll (coarseCoeff, i)
        {
            activeCoarseCoeff[i] = pTraits<scalarType>::zero;
        }

        forAll (Coeff, i)
        {
            activeCoarseCoeff[agglomIndex_[i]] += activeCoeff[i];
        }
    }
    else
    {
        FatalErrorIn
        (
            "void BlockMatrixAgglomeration<Type>::restrictDiagDecoupled()"
            " const"
        )   << "Problem in coeff type morphing"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockMatrixAgglomeration<Type>::BlockMatrixAgglomeration
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
    agglomIndex_(matrix_.lduAddr().size()),
    groupSize_(groupSize),
    nSolo_(0),
    nCoarseEqns_(0),
    coarsen_(false),
    lTime_()
{
    calcAgglomeration();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockMatrixAgglomeration<Type>::~BlockMatrixAgglomeration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockAMGLevel<Type> >
Foam::BlockMatrixAgglomeration<Type>::restrictMatrix() const
{
    if (!coarsen_)
    {
        FatalErrorIn
        (
            "autoPtr<BlockAMGLevel<Type> > "
            "BlockMatrixAgglomeration<Type>::restrictMatrix() const"
        )   << "Requesting coarse matrix when it cannot be created"
            << abort(FatalError);
    }

    // Construct the coarse matrix and ldu addressing for the next level
    // Algorithm:
    // 1) Loop through all fine coeffs.  If the agglom labels on two sides are
    //    different, this creates a coarse coeff. Define owner and neighbour
    //    for this coeff based on cluster IDs.
    // 2) Check if the coeff has been seen before. If yes, add the coefficient
    //    to the appropriate field (stored with the equation). If no, create
    //    a new coeff with neighbour ID and add the coefficient
    // 3) Once all the coeffs have been created, loop through all clusters and
    //    insert the coeffs in the upper order. At the same time, collect the
    //    owner and neighbour addressing.
    // 4) Agglomerate the diagonal by summing up the fine diagonal

    // Get addressing
    const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();

    const label nFineCoeffs = upperAddr.size();

#   ifdef FULLDEBUG
    if (agglomIndex_.size() != matrix_.lduAddr().size())
    {
        FatalErrorIn
        (
            "autoPtr<BlockLduMatrix<Type> >"
            "BlockMatrixAgglomeration<Type>::restrictMatrix() const"
        )   << "agglomIndex array does not correspond to fine level. " << endl
            << " Size: " << agglomIndex_.size()
            << " number of equations: " << matrix_.lduAddr().size()
            << abort(FatalError);
    }
#   endif

    // Does the matrix have solo equations
    bool soloEqns = nSolo_ > 0;

    // Storage for block neighbours and coefficients

    // Guess initial maximum number of neighbours in block
    label maxNnbrs = 10;

    // Number of neighbours per block
    labelList blockNnbrs(nCoarseEqns_, 0);

    // Setup initial packed storage for neighbours and coefficients
    labelList blockNbrsData(maxNnbrs*nCoarseEqns_);

    // Create face-restriction addressing
    labelList coeffRestrictAddr(nFineCoeffs);

    // Initial neighbour array (not in upper-triangle order)
    labelList initCoarseNeighb(nFineCoeffs);

    // Counter for coarse coeffs
    label nCoarseCoeffs = 0;

    // Loop through all fine coeffs
    forAll (upperAddr, fineCoeffi)
    {
        label rmUpperAddr = agglomIndex_[upperAddr[fineCoeffi]];
        label rmLowerAddr = agglomIndex_[lowerAddr[fineCoeffi]];

        // If the coefficient touches block zero and solo equations are
        // present, skip it
        if (soloEqns && (rmUpperAddr == 0 || rmLowerAddr == 0))
        {
            continue;
        }

        if (rmUpperAddr == rmLowerAddr)
        {
            // For each fine coeff inside of a coarse cluster keep the address
            // of the cluster corresponding to the coeff in the
            // coeffRestrictAddr as a negative index
            coeffRestrictAddr[fineCoeffi] = -(rmUpperAddr + 1);
        }
        else
        {
            // This coeff is a part of a coarse coeff

            label cOwn = rmUpperAddr;
            label cNei = rmLowerAddr;

            // Get coarse owner and neighbour
            if (rmUpperAddr > rmLowerAddr)
            {
                cOwn = rmLowerAddr;
                cNei = rmUpperAddr;
            }

            // Check the neighbour to see if this coeff has already been found
            bool nbrFound = false;
            label& ccnCoeffs = blockNnbrs[cOwn];

            for (int i = 0; i < ccnCoeffs; i++)
            {
                if (initCoarseNeighb[blockNbrsData[maxNnbrs*cOwn + i]] == cNei)
                {
                    nbrFound = true;
                    coeffRestrictAddr[fineCoeffi] =
                        blockNbrsData[maxNnbrs*cOwn + i];
                    break;
                }
            }

            if (!nbrFound)
            {
                if (ccnCoeffs >= maxNnbrs)
                {
                    // Double the size of list and copy data
                    label oldMaxNnbrs = maxNnbrs;
                    maxNnbrs *= 2;

                    // Resize and copy list
                    const labelList oldBlockNbrsData = blockNbrsData;
                    blockNbrsData.setSize(maxNnbrs*nCoarseEqns_);

                    forAll (blockNnbrs, i)
                    {
                        for (int j = 0; j < blockNnbrs[i]; j++)
                        {
                            blockNbrsData[maxNnbrs*i + j] =
                                oldBlockNbrsData[oldMaxNnbrs*i + j];
                        }
                    }
                }

                blockNbrsData[maxNnbrs*cOwn + ccnCoeffs] = nCoarseCoeffs;
                initCoarseNeighb[nCoarseCoeffs] = cNei;
                coeffRestrictAddr[fineCoeffi] = nCoarseCoeffs;
                ccnCoeffs++;

                // New coarse coeff created
                nCoarseCoeffs++;
            }
        }
    } // End for all fine coeffs


    // Renumber into upper-triangular order

    // All coarse owner-neighbour storage
    labelList coarseOwner(nCoarseCoeffs);
    labelList coarseNeighbour(nCoarseCoeffs);
    labelList coarseCoeffMap(nCoarseCoeffs);

    label coarseCoeffi = 0;

    forAll (blockNnbrs, cci)
    {
        label* cCoeffs = &blockNbrsData[maxNnbrs*cci];
        label ccnCoeffs = blockNnbrs[cci];

        for (int i = 0; i < ccnCoeffs; i++)
        {
            coarseOwner[coarseCoeffi] = cci;
            coarseNeighbour[coarseCoeffi] = initCoarseNeighb[cCoeffs[i]];
            coarseCoeffMap[cCoeffs[i]] = coarseCoeffi;
            coarseCoeffi++;
        }
    }

    forAll(coeffRestrictAddr, fineCoeffi)
    {
        label rmUpperAddr = agglomIndex_[upperAddr[fineCoeffi]];
        label rmLowerAddr = agglomIndex_[lowerAddr[fineCoeffi]];

        // If the coefficient touches block zero and solo equations are
        // present, skip it
        if (soloEqns && (rmUpperAddr == 0 || rmLowerAddr == 0))
        {
            continue;
        }

        if (coeffRestrictAddr[fineCoeffi] >= 0)
        {
            coeffRestrictAddr[fineCoeffi] =
                coarseCoeffMap[coeffRestrictAddr[fineCoeffi]];
        }
    }

    // Clear the temporary storage for the coarse matrix data
    blockNnbrs.setSize(0);
    blockNbrsData.setSize(0);
    initCoarseNeighb.setSize(0);
    coarseCoeffMap.setSize(0);


    // Create coarse-level coupled interfaces

    // Create coarse interfaces, addressing and coefficients
    const label interfaceSize =
        const_cast<BlockLduMatrix<Type>& >(matrix_).interfaces().size();

    const typename BlockLduInterfaceFieldPtrsList<Type>::Type&
        interfaceFields =
        const_cast<BlockLduMatrix<Type>&>(matrix_).interfaces();

    // Set the coarse interfaces and coefficients
    lduInterfacePtrsList coarseInterfaces(interfaceSize);

    labelListList coarseInterfaceAddr(interfaceSize);

    // Add the coarse level

    // Set the coarse ldu addressing onto the list
    autoPtr<lduPrimitiveMesh> coarseAddrPtr
    (
        new lduPrimitiveMesh
        (
            nCoarseEqns_,
            coarseOwner,
            coarseNeighbour,
            Pstream::worldComm, //HJ, AMG Comm fineMesh.comm(),
            true
        )
    );

    // Initialise transfer of restrict addressing on the interface
    // HJ, consider blocking comms.  HJ, 9/Jun/2016
    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            interfaceFields[intI].coupledInterface().initInternalFieldTransfer
            (
                Pstream::blocking,
                agglomIndex_
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
                        agglomIndex_
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
                    fineInterface.interfaceInternalField(agglomIndex_),
                    fineInterfaceAddr[intI]
                ).ptr()
            );
        }
    }

    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            const AMGInterface& coarseInterface =
                refCast<const AMGInterface>(coarseInterfaces[intI]);

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

    // Set the coarse level matrix
    autoPtr<BlockLduMatrix<Type> > coarseMatrixPtr
    (
        new BlockLduMatrix<Type>(coarseAddrPtr())
    );
    BlockLduMatrix<Type>& coarseMatrix = coarseMatrixPtr();

    // Get interfaces from coarse matrix
    typename BlockLduInterfaceFieldPtrsList<Type>::Type&
        coarseInterfaceFieldsTransfer = coarseMatrix.interfaces();

    // Aggolmerate the upper and lower coupled coefficients
    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            const AMGInterface& coarseInterface =
                refCast<const AMGInterface>(coarseInterfaces[intI]);

            coarseInterfaceFieldsTransfer.set
            (
                intI,
                BlockAMGInterfaceField<Type>::New
                (
                    coarseInterface,
                    interfaceFields[intI]
                ).ptr()
            );

            // Since the type of agglomeration is now templated, agglomeration
            // of block coefficients must be done by a FIELD (not interface)
            // via a new set of virtual functions
            // HJ, 16/Mar/2016

            // Note: in the scalar AMG, agglomeration is done by the interface
            // (always scalar) but in the block matrix it is done by a
            // templated block interface field
            // HJ, 16/Mar/2016

            // Cast the interface into AMG type
            const BlockAMGInterfaceField<Type>& coarseField =
                refCast<const BlockAMGInterfaceField<Type> >
                (
                    coarseInterfaceFieldsTransfer[intI]
                );

            coarseMatrix.coupleUpper().set
            (
                intI,
                coarseField.agglomerateBlockCoeffs
                (
                    matrix_.coupleUpper()[intI]
                )
            );

            coarseMatrix.coupleLower().set
            (
                intI,
                coarseField.agglomerateBlockCoeffs
                (
                    matrix_.coupleLower()[intI]
                )
            );
        }
    }

    // Matrix restriction done!

    typedef CoeffField<Type> TypeCoeffField;

    typedef typename TypeCoeffField::squareTypeField squareTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;

    TypeCoeffField& coarseUpper = coarseMatrix.upper();
    TypeCoeffField& coarseDiag = coarseMatrix.diag();
    const TypeCoeffField& fineUpper = matrix_.upper();
    const TypeCoeffField& fineDiag = matrix_.diag();

    // KRJ: 2013-01-31: Many cases needed as there are different combinations

    // Note:
    // In coarsening, the off-diagonal coefficient type should be preserved
    // and the case of taking out of off-diag from diag may need to be handled
    // separately (expand the off-diag coeff to diag type before removing it
    // from the diag coefficient).  Since this has not been encountered yet
    // only matching diag/off-diag types are handled.
    // HJ, 15/Feb/2016
    if (matrix_.symmetric())
    {
        if
        (
            fineDiag.activeType() == blockCoeffBase::SQUARE
         || fineUpper.activeType() == blockCoeffBase::SQUARE
        )
        {
            squareTypeField& activeCoarseDiag = coarseDiag.asSquare();

            squareTypeField& activeCoarseUpper = coarseUpper.asSquare();
            const squareTypeField& activeFineUpper = fineUpper.asSquare();

            // Use lower as transpose of upper
            squareTypeField activeFineUpperTranspose =
                activeFineUpper.T();

            restrictDiag(fineDiag, coarseDiag);

            agglomerateCoeffs
            (
                coeffRestrictAddr,
                activeCoarseDiag,
                activeCoarseUpper,
                activeFineUpper,
                activeFineUpperTranspose
            );
        }
        else if
        (
            fineDiag.activeType() == blockCoeffBase::LINEAR
         || fineUpper.activeType() == blockCoeffBase::LINEAR
        )
        {
            linearTypeField& activeCoarseDiag = coarseDiag.asLinear();

            linearTypeField& activeCoarseUpper = coarseUpper.asLinear();
            const linearTypeField& activeFineUpper = fineUpper.asLinear();

            restrictDiag(fineDiag, coarseDiag);

            agglomerateCoeffs
            (
                coeffRestrictAddr,
                activeCoarseDiag,
                activeCoarseUpper,
                activeFineUpper,
                activeFineUpper
            );
        }
        else if
        (
            fineDiag.activeType() == blockCoeffBase::SCALAR
         || fineUpper.activeType() == blockCoeffBase::SCALAR
        )
        {
            scalarTypeField& activeCoarseDiag = coarseDiag.asScalar();

            scalarTypeField& activeCoarseUpper = coarseUpper.asScalar();
            const scalarTypeField& activeFineUpper = fineUpper.asScalar();

            restrictDiag(fineDiag, coarseDiag);

            agglomerateCoeffs
            (
                coeffRestrictAddr,
                activeCoarseDiag,
                activeCoarseUpper,
                activeFineUpper,
                activeFineUpper
            );
        }
        else
        {
            FatalErrorIn
            (
                "autoPtr<BlockAMGLevel<Type> >"
                "BlockMatrixAgglomeration<Type>::restrictMatrix() const"
            )   << "Matrix coeff type morphing error, symmetric matrix"
                << abort(FatalError);
        }
    }
    else // asymmetric matrix
    {
        TypeCoeffField& coarseLower = coarseMatrix.lower();
        const TypeCoeffField& fineLower = matrix_.lower();

        if
        (
            fineDiag.activeType() == blockCoeffBase::SQUARE
         || fineUpper.activeType() == blockCoeffBase::SQUARE
        )
        {
            squareTypeField& activeCoarseDiag = coarseDiag.asSquare();

            squareTypeField& activeCoarseUpper = coarseUpper.asSquare();
            const squareTypeField& activeFineUpper = fineUpper.asSquare();

            squareTypeField& activeCoarseLower = coarseLower.asSquare();
            const squareTypeField& activeFineLower = fineLower.asSquare();

            restrictDiag(fineDiag, coarseDiag);

            agglomerateCoeffs
            (
                coeffRestrictAddr,
                activeCoarseDiag,
                activeCoarseUpper,
                activeFineUpper,
                activeCoarseLower,
                activeFineLower
            );
        }
        else if
        (
            fineDiag.activeType() == blockCoeffBase::LINEAR
         || fineUpper.activeType() == blockCoeffBase::LINEAR
        )
        {
            linearTypeField& activeCoarseDiag = coarseDiag.asLinear();

            linearTypeField& activeCoarseUpper = coarseUpper.asLinear();
            const linearTypeField& activeFineUpper = fineUpper.asLinear();

            linearTypeField& activeCoarseLower = coarseLower.asLinear();
            const linearTypeField& activeFineLower = fineLower.asLinear();

            restrictDiag(fineDiag, coarseDiag);

            agglomerateCoeffs
            (
                coeffRestrictAddr,
                activeCoarseDiag,
                activeCoarseUpper,
                activeFineUpper,
                activeCoarseLower,
                activeFineLower
            );
        }
        else if
        (
            fineDiag.activeType() == blockCoeffBase::SCALAR
         || fineUpper.activeType() == blockCoeffBase::SCALAR
        )
        {
            scalarTypeField& activeCoarseDiag = coarseDiag.asScalar();

            scalarTypeField& activeCoarseUpper = coarseUpper.asScalar();
            const scalarTypeField& activeFineUpper = fineUpper.asScalar();

            scalarTypeField& activeCoarseLower = coarseLower.asScalar();
            const scalarTypeField& activeFineLower = fineLower.asScalar();

            restrictDiag(fineDiag, coarseDiag);

            agglomerateCoeffs
            (
                coeffRestrictAddr,
                activeCoarseDiag,
                activeCoarseUpper,
                activeFineUpper,
                activeCoarseLower,
                activeFineLower
            );
        }
        else
        {
            FatalErrorIn
            (
                "autoPtr<BlockAMGLevel<Type> >"
                "BlockMatrixAgglomeration<Type>::restrictMatrix() const"
            )   << "Matrix coeff type morphing error, asymmetric matrix"
                << abort(FatalError);
        }
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
            this->groupSize(),
            this->minCoarseEqns()
        )
    );
}


template<class Type>
void Foam::BlockMatrixAgglomeration<Type>::restrictResidual
(
    const Field<Type>& res,
    Field<Type>& coarseRes
) const
{
    coarseRes = pTraits<Type>::zero;

    forAll (res, i)
    {
        coarseRes[agglomIndex_[i]] += res[i];
    }
}


template<class Type>
void Foam::BlockMatrixAgglomeration<Type>::prolongateCorrection
(
    Field<Type>& x,
    const Field<Type>& coarseX
) const
{
    forAll (x, i)
    {
        x[i] += coarseX[agglomIndex_[i]];
    }
}


// ************************************************************************* //
