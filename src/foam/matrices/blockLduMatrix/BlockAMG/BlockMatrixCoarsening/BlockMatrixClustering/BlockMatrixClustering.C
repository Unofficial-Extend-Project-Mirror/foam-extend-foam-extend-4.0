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
    BlockMatrixClustering

Description
    Block matrix AMG coarsening by Jasak clustering algorithm

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "BlockMatrixClustering.H"
#include "boolList.H"
#include "tolerancesSwitch.H"
#include "coeffFields.H"
#include "BlockAMGInterfaceField.H"
#include "coarseBlockAMGLevel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::debug::tolerancesSwitch
Foam::BlockMatrixClustering<Type>::weightFactor_
(
    "aamgWeightFactor",
    0.65
);


template<class Type>
const Foam::debug::tolerancesSwitch
Foam::BlockMatrixClustering<Type>::diagFactor_
(
    "aamgDiagFactor",
    1e-8
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::BlockMatrixClustering<Type>::calcClustering()
{
    if (matrix_.diagonal())
    {
        // Diag only matrix.  Reset and return
        agglomIndex_ = 0;
        nCoarseEqns_ = 1;

        return;
    }

    // Algorithm:
    // 1) Calculate appropriate connection strength for symmetric and asymmetric
    //    matrix
    // 2) Collect solo (disconnected) cells and remove them from agglomeration
    //    by placing them into cluster zero (solo group)
    // 3) Loop through all equations and for each equation find the best fit
    //    neighbour.  If all neighbours are grouped, add equation
    //    to best group

    // Initialise child array
    agglomIndex_ = -1;

    const label nRows = matrix_.lduAddr().size();

    // Get matrix addressing
    const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();
    const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
    const unallocLabelList& losortAddr = matrix_.lduAddr().losortAddr();

    const unallocLabelList& ownerStartAddr =
        matrix_.lduAddr().ownerStartAddr();
    const unallocLabelList& losortStartAddr =
        matrix_.lduAddr().losortStartAddr();


    // Calculate clustering

    // Get matrix coefficients and norms  Note: the norm
    // may be signed, ie. it will take the sign of the coefficient

    const CoeffField<Type>& diag = matrix_.diag();
    scalarField normDiag(diag.size());
    normPtr_->normalize(normDiag, diag);

    scalarField normUpper(upperAddr.size(), 0);
    scalarField normLower(upperAddr.size(), 0);

    // Note:
    // Matrix properties are no longer assumed, eg. the diag and off-diag
    // sign is checked
    // HJ, 29/Mar/2017

    // Note: negative connections are eliminated in max(...) below
    // HJ, 30/Mar/2017

    if (matrix_.thereIsUpper())
    {
        normPtr_->normalize(normUpper, matrix_.upper());

        // Owner: upper triangle
        forAll (lowerAddr, coeffI)
        {
            // Sign of strong positive upper is opposite of the sign of
            // its diagonal coefficient
            normUpper[coeffI] = Foam::max
            (
                -1*sign(normDiag[lowerAddr[coeffI]])*normUpper[coeffI],
                0
            );
        }
    }

    if (matrix_.thereIsLower())
    {
        normPtr_->normalize(normLower, matrix_.lower());

        // Neighbour: lower triangle
        forAll (lowerAddr, coeffI)
        {
            // Sign of strong positive upper is opposite of the sign of
            // its diagonal coefficient
            normLower[coeffI] = Foam::max
            (
                -1*sign(normDiag[upperAddr[coeffI]])*normLower[coeffI],
                0
            );
        }
    }
    else
    {
        normLower = normUpper;
    }

    // Take magnitude of the diagonal norm after the normalisation
    // of upper and lower is complete
    // Note: upper and lower are already normalised
    normDiag = mag(normDiag);

    labelList sizeOfGroups(nRows, 0);

    nCoarseEqns_ = 0;

    // Gather disconnected and weakly connected equations into cluster zero
    // Weak connection is assumed to be the one where the off-diagonal
    // coefficient is smaller than diagFactor_()*diag
    {
        // Algorithm
        // Mark all cells to belong to zero cluster
        // Go through all upper and lower coefficients and for the ones
        // larger than threshold mark the equations out of cluster zero

        scalarField scaledNormDiag = diagFactor_()*normDiag;

        boolList zeroCluster(normDiag.size(), true);

        if (matrix_.symmetric())
        {
            // Owner: upper triangle
            forAll (lowerAddr, coeffI)
            {
                if (normUpper[coeffI] > scaledNormDiag[lowerAddr[coeffI]])
                {
                    zeroCluster[lowerAddr[coeffI]] = false;
                }
            }

            // Neighbour: lower triangle with symm coefficients
            forAll (upperAddr, coeffI)
            {
                if (normUpper[coeffI] > scaledNormDiag[upperAddr[coeffI]])
                {
                    zeroCluster[upperAddr[coeffI]] = false;
                }
            }
        }
        else if (matrix_.asymmetric())
        {
            // Owner: upper triangle
            forAll (lowerAddr, coeffI)
            {
                if (normUpper[coeffI] > scaledNormDiag[lowerAddr[coeffI]])
                {
                    zeroCluster[lowerAddr[coeffI]] = false;
                }
            }

            // Neighbour: lower triangle with lower coeffs
            forAll (upperAddr, coeffI)
            {
                if (normLower[coeffI] > scaledNormDiag[upperAddr[coeffI]])
                {
                    zeroCluster[upperAddr[coeffI]] = false;
                }
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
            sizeOfGroups[nCoarseEqns_] = nSolo_;

            nCoarseEqns_++;
        }
    }


    // Loop through cells
    // - if the cell is not already grouped, open a new group (seed)
    // - find stroungest grouped and ungrouped connection:
    //   - grouped connection is towards an already grouped cell
    //   - ungrouped connection is towards an ungrouped cell
    // - decide if a grouped or ungrouped connection is better

    label curEqn, indexUngrouped, indexGrouped, colI, groupPassI,
        nextUngrouped, nextGrouped;

    scalar curWeight, weightUngrouped, weightGrouped;

    for (label rowI = 0; rowI < nRows; rowI++)
    {
        if (agglomIndex_[rowI] == -1)
        {
            // Found new ungrouped equation
            curEqn = rowI;

            // Reset grouped and upgrouped index
            indexUngrouped = -1;
            indexGrouped = -1;

            // Make next group (coarse equation) and search for neighbours
            agglomIndex_[curEqn] = nCoarseEqns_;

            // Work on the group until the min group size is satisfied
            // As each new element of the group is found, group search starts
            // from this element until group cluster size is satisfied
            for
            (
                groupPassI = 1;
                groupPassI < minGroupSize_;
                groupPassI++
            )
            {
                weightUngrouped = 0;
                weightGrouped = 0;

                indexUngrouped = -1;
                indexGrouped = -1;

                nextUngrouped = -1;
                nextGrouped = -1;

                // Visit all neighbour equations

                // Visit upper triangle coefficients from the equation
                for
                (
                    label rowCoeffI = ownerStartAddr[curEqn];
                    rowCoeffI < ownerStartAddr[curEqn + 1];
                    rowCoeffI++
                )
                {
                    // Get column index, upper triangle
                    colI = upperAddr[rowCoeffI];

                    // curWeight = normUpper[rowCoeffI]/
                    //     // normDiag[curEqn];
                    //     max(normDiag[curEqn], normDiag[colI]);

                    curWeight = Foam::min
                    (
                        normUpper[rowCoeffI]/normDiag[curEqn],
                        normLower[rowCoeffI]/normDiag[colI]
                    );

                    if (agglomIndex_[colI] == -1)
                    {
                        // Found ungrouped neighbour
                        if
                        (
                            indexUngrouped == -1
                         || curWeight > weightUngrouped
                        )
                        {
                            // Found or updated ungrouped neighbour
                            indexUngrouped = rowCoeffI;
                            weightUngrouped = curWeight;

                            // Record possible next equation for search
                            nextUngrouped = colI;
                        }
                    }
                    else if (agglomIndex_[curEqn] != agglomIndex_[colI])
                    {
                        // Check for neighbour in solo group
                        if (nSolo_ == 0 || agglomIndex_[colI] != 0)
                        {
                            // Found neighbour belonging to other group
                            if (indexGrouped == -1 || curWeight > weightGrouped)
                            {
                                // Found or updated grouped neighbour
                                indexGrouped = rowCoeffI;
                                weightGrouped = curWeight;

                                // Record possible next equation for search
                                nextGrouped = colI;
                            }
                        }
                    }
                }

                // Visit lower triangle coefficients from the equation
                for
                (
                    label rowCoeffI = losortStartAddr[curEqn];
                    rowCoeffI < losortStartAddr[curEqn + 1];
                    rowCoeffI++
                )
                {
                    // Get column index, lower triangle
                    colI = lowerAddr[losortAddr[rowCoeffI]];

                    // curWeight = normLower[losortAddr[rowCoeffI]]/
                    //     // normDiag[curEqn];
                    //     max(normDiag[curEqn], normDiag[colI]);

                    curWeight = Foam::min
                    (
                        normLower[losortAddr[rowCoeffI]]/normDiag[curEqn],
                        normUpper[losortAddr[rowCoeffI]]/normDiag[colI]
                    );

                    if (agglomIndex_[colI] == -1)
                    {
                        // Found first or better ungrouped neighbour
                        if
                        (
                            indexUngrouped == -1
                         || curWeight > weightUngrouped
                        )
                        {
                            // Found or updated ungrouped neighbour
                            indexUngrouped = rowCoeffI;
                            weightUngrouped = curWeight;

                            // Record possible next equation for search
                            nextUngrouped = colI;
                        }
                    }
                    else if (agglomIndex_[curEqn] != agglomIndex_[colI])
                    {
                        // Check for neighbour in solo group
                        if (nSolo_ == 0 || agglomIndex_[colI] != 0)
                        {
                            // Found first of better neighbour belonging to
                            // other group
                            if
                            (
                                indexGrouped == -1
                             || curWeight > weightGrouped
                            )
                            {
                                // Found or updated grouped neighbour
                                indexGrouped = rowCoeffI;
                                weightGrouped = curWeight;

                                // Record possible next equation for search
                                nextGrouped = colI;
                            }
                        }
                    }
                }

                // Decide what to do with current equation

                // If ungrouped candidate exists and it is stronger
                // than best weighted grouped candidate, use it
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

                    agglomIndex_[nextUngrouped] = agglomIndex_[curEqn];
                    sizeOfGroups[agglomIndex_[curEqn]]++;

                    // Search from nextEqn
                    curEqn = nextUngrouped;
                }
                else
                {
                    // Group full or cannot be extended with new
                    // ungrouped candidates
                    break;
                }
            }


            // Finished group passes
            if
            (
                groupPassI > 1
             || indexGrouped == -1
             || sizeOfGroups[agglomIndex_[nextGrouped]] >= maxGroupSize_
            )
            {
                // If this is a solo cell and a group is available
                // force it into the group irrespective of size
                if
                (
                    sizeOfGroups[agglomIndex_[rowI]] == 0
                 && indexGrouped != -1

                )
                {
                    // Group exists, but it's too big.  Add it anyway
                    agglomIndex_[rowI] = agglomIndex_[nextGrouped];
                    sizeOfGroups[agglomIndex_[nextGrouped]]++;
                }
                else
                {
                    // No group and no solo group.  Make its own group
                    sizeOfGroups[agglomIndex_[rowI]]++;
                    nCoarseEqns_++;
                }
            }
            else
            {
                // Dump current cell into the best group available
                agglomIndex_[rowI] = agglomIndex_[nextGrouped];
                sizeOfGroups[agglomIndex_[nextGrouped]]++;
            }
        }
    }

    // Resize group size count
    sizeOfGroups.setSize(nCoarseEqns_);

    // The decision on parallel agglomeration needs to be made for the
    // whole gang of processes; otherwise I may end up with a different
    // number of agglomeration levels on different processors.

    // If the number of coarse equations is less than minimum and
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

    if (blockLduMatrix::debug >= 3)
    {
        // Count singleton clusters
        label nSingleClusters = 0;

        // Adjust start based on solo cells cluster status
        label start = 0;

        if (nSolo_ > 0)
        {
            start = 1;
        }

        for (label i = start; i < sizeOfGroups.size(); i++)
        {
            if (sizeOfGroups[i] == 1)
            {
                nSingleClusters++;
            }
        }

        Pout<< "Coarse level size: " << nCoarseEqns_
            << " nSolo = " << nSolo_
            << " cluster (" << min(sizeOfGroups) << " "
            << max(sizeOfGroups)
            << ").  N singleton clusters = " << nSingleClusters;

        if (coarsen_)
        {
            Pout << ".  Accepted" << endl;
        }
        else
        {
            Pout << ".  Rejected" << endl;
        }
    }
}


template<class Type>
void Foam::BlockMatrixClustering<Type>::restrictDiag
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

        // Reset coefficients to zero
        activeCoarseCoeff = pTraits<squareType>::zero;

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

        // Reset coefficients to zero
        activeCoarseCoeff = pTraits<linearType>::zero;

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

        // Reset coefficients to zero
        activeCoarseCoeff = pTraits<scalarType>::zero;

        forAll (Coeff, i)
        {
            activeCoarseCoeff[agglomIndex_[i]] += activeCoeff[i];
        }
    }
    else
    {
        FatalErrorIn
        (
            "void  BlockMatrixClustering<Type>::restrictDiag() const"
        )   << "Problem in coeff type morphing"
            << abort(FatalError);
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockMatrixClustering<Type>::agglomerateCoeffs
(
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

    // Reset coefficients to zero.  Cannot touch the diagonal
    activeCoarseUpper = pTraits<ULType>::zero;

    forAll(coeffRestrictAddr_, fineCoeffI)
    {
        label rmUpperAddr = agglomIndex_[upperAddr[fineCoeffI]];
        label rmLowerAddr = agglomIndex_[lowerAddr[fineCoeffI]];

        // If the coefficient touches block zero and
        //  solo equations are present, skip it
        if (soloEqns && (rmUpperAddr == 0 || rmLowerAddr == 0))
        {
            continue;
        }

        label cCoeff = coeffRestrictAddr_[fineCoeffI];

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
void Foam::BlockMatrixClustering<Type>::agglomerateCoeffs
(
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

    // Reset coefficients to zero.  Cannot touch the diagonal
    activeCoarseUpper = pTraits<ULType>::zero;
    activeCoarseLower = pTraits<ULType>::zero;

    forAll(coeffRestrictAddr_, fineCoeffI)
    {
        label rmUpperAddr = agglomIndex_[upperAddr[fineCoeffI]];
        label rmLowerAddr = agglomIndex_[lowerAddr[fineCoeffI]];

        // If the coefficient touches block zero and
        //  solo equations are present, skip it
        if (soloEqns && (rmUpperAddr == 0 || rmLowerAddr == 0))
        {
            continue;
        }

        label cCoeff = coeffRestrictAddr_[fineCoeffI];

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
void Foam::BlockMatrixClustering<Type>::restrictDiagDecoupled
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

        // Reset coefficients to zero
        activeCoarseCoeff = pTraits<linearType>::zero;

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

        // Reset coefficients to zero
        activeCoarseCoeff = pTraits<scalarType>::zero;

        forAll (Coeff, i)
        {
            activeCoarseCoeff[agglomIndex_[i]] += activeCoeff[i];
        }
    }
    else
    {
        FatalErrorIn
        (
            "void BlockMatrixClustering<Type>::restrictDiagDecoupled()"
            " const"
        )   << "Problem in coeff type morphing"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockMatrixClustering<Type>::BlockMatrixClustering
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict,
    const label groupSize,
    const label minCoarseEqns
)
:
    BlockMatrixCoarsening<Type>(matrix, dict, groupSize, minCoarseEqns),
    matrix_(matrix),
    minGroupSize_(readLabel(dict.lookup("minGroupSize"))),
    maxGroupSize_(readLabel(dict.lookup("maxGroupSize"))),
    normPtr_(BlockCoeffNorm<Type>::New(dict)),
    agglomIndex_(matrix_.lduAddr().size()),
    coeffRestrictAddr_(),
    groupSize_(groupSize),
    nSolo_(0),
    nCoarseEqns_(0),
    coarsen_(false)
{
    calcClustering();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockMatrixClustering<Type>::~BlockMatrixClustering()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockAMGLevel<Type> >
Foam::BlockMatrixClustering<Type>::restrictMatrix() const
{
    if (!coarsen_)
    {
        FatalErrorIn
        (
            "autoPtr<BlockAMGLevel<Type> > "
            "BlockMatrixClustering<Type>::restrictMatrix() const"
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

#   ifdef FULLDEBUG
    if (agglomIndex_.size() != matrix_.lduAddr().size())
    {
        FatalErrorIn
        (
            "autoPtr<BlockLduMatrix<Type> >"
            "BlockMatrixClustering<Type>::restrictMatrix() const"
        )   << "agglomIndex array does not correspond to fine level. " << endl
            << " Size: " << agglomIndex_.size()
            << " number of equations: " << matrix_.lduAddr().size()
            << abort(FatalError);
    }
#   endif

    // If the matrix will be coarsened, create off-diagonal agglomeration
        // Does the matrix have solo equations
    bool soloEqns = nSolo_ > 0;

    // Storage for block neighbours and coefficients

    // Guess initial maximum number of neighbours in block
    label maxNnbrs = 10;

    // Number of neighbours per block
    labelList blockNnbrs(nCoarseEqns_, 0);

    // Setup initial packed storage for neighbours and coefficients
    labelList blockNbrsData(maxNnbrs*nCoarseEqns_);

    const label nFineCoeffs = upperAddr.size();

    // Create face-restriction addressing
    coeffRestrictAddr_.setSize(nFineCoeffs);

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
            // coeffRestrictAddr_ as a negative index
            coeffRestrictAddr_[fineCoeffi] = -(rmUpperAddr + 1);
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
                    coeffRestrictAddr_[fineCoeffi] =
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
                coeffRestrictAddr_[fineCoeffi] = nCoarseCoeffs;
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

    forAll(coeffRestrictAddr_, fineCoeffi)
    {
        label rmUpperAddr = agglomIndex_[upperAddr[fineCoeffi]];
        label rmLowerAddr = agglomIndex_[lowerAddr[fineCoeffi]];

        // If the coefficient touches block zero and solo equations are
        // present, skip it
        if (soloEqns && (rmUpperAddr == 0 || rmLowerAddr == 0))
        {
            continue;
        }

        if (coeffRestrictAddr_[fineCoeffi] >= 0)
        {
            coeffRestrictAddr_[fineCoeffi] =
                coarseCoeffMap[coeffRestrictAddr_[fineCoeffi]];
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
        interfaceFields = matrix_.interfaces();

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
            Pstream::worldComm,
            true
        )
    );

    // Initialise transfer of restrict addressing on the interface
    // HJ, must use blocking comms.  HJ, 9/Jun/2016
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
            this->groupSize(),
            this->minCoarseEqns()
        )
    );
}


template<class Type>
void Foam::BlockMatrixClustering<Type>::updateMatrix
(
    BlockLduMatrix<Type>& coarseMatrix
) const
{
    // Get interfaces from fine matrix
    const typename BlockLduInterfaceFieldPtrsList<Type>::Type&
        interfaceFields = matrix_.interfaces();

    // Get interfaces from coarse matrix
    lduInterfacePtrsList coarseInterfaces = coarseMatrix.mesh().interfaces();

    // Get interfaces fields from coarse matrix
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

            // Since the type of clustering is now templated, clustering
            // of block coefficients must be done by a FIELD (not interface)
            // via a new set of virtual functions
            // HJ, 16/Mar/2016

            // Note: in the scalar AMG, clustering is done by the interface
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
                "BlockMatrixClustering<Type>::restrictMatrix() const"
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
                "BlockMatrixClustering<Type>::restrictMatrix() const"
            )   << "Matrix coeff type morphing error, asymmetric matrix"
                << abort(FatalError);
        }
    }
}


template<class Type>
void Foam::BlockMatrixClustering<Type>::restrictResidual
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
void Foam::BlockMatrixClustering<Type>::prolongateCorrection
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
