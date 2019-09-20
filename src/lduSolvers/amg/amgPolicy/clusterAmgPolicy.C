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
    clusterAmgPolicy

Description
    Block matrix AMG coarsening by Jasak clustering algorithm

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "clusterAmgPolicy.H"
#include "amgMatrix.H"
#include "boolList.H"
#include "addToRunTimeSelectionTable.H"
#include "AMGInterfaceField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(clusterAmgPolicy, 0);

    addToRunTimeSelectionTable(amgPolicy, clusterAmgPolicy, matrix);

} // End namespace Foam


const Foam::debug::tolerancesSwitch Foam::clusterAmgPolicy::weightFactor_
(
    "aamgWeightFactor",
    0.65
);


const Foam::debug::tolerancesSwitch Foam::clusterAmgPolicy::diagFactor_
(
    "aamgDiagFactor",
    1e-8
);


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

void Foam::clusterAmgPolicy::calcChild()
{
    if (matrix().diagonal())
    {
        // Diag only matrix.  Reset and return
        child_ = 0;
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
    child_ = -1;

    const label nRows = matrix().lduAddr().size();

    // Get matrix addressing
    const unallocLabelList& lowerAddr = matrix().lduAddr().lowerAddr();
    const unallocLabelList& upperAddr = matrix().lduAddr().upperAddr();
    const unallocLabelList& losortAddr = matrix().lduAddr().losortAddr();

    const unallocLabelList& ownerStartAddr =
        matrix().lduAddr().ownerStartAddr();
    const unallocLabelList& losortStartAddr =
        matrix().lduAddr().losortStartAddr();


    // Calculate agglomeration

    // Get magnitudes of matrix coefficients

    const scalarField& diag = matrix().diag();
    const scalarField magDiag = mag(matrix().diag());

    // Calculate magnitude of strong positive connections
    scalarField magUpper(upperAddr.size(), 0);
    scalarField magLower(upperAddr.size(), 0);

    // Note:
    // Matrix properties are no longer assumed, eg. the diag and off-diag
    // sign is checked
    // HJ, 29/Mar/2017

    // Note: negative connections are eliminated in max(...) below
    // HJ, 30/Mar/2017

    if (matrix().hasUpper())
    {
        const scalarField& upper = matrix().upper();

        // Owner: upper triangle
        forAll (lowerAddr, coeffI)
        {
            // Sign of strong positive upper is opposite of the sign of
            // its diagonal coefficient
            magUpper[coeffI] =
                Foam::max(-1*sign(diag[lowerAddr[coeffI]])*upper[coeffI], 0);
        }
    }

    if (matrix().hasLower())
    {
        const scalarField& lower = matrix().lower();

        // Neighbour: lower triangle
        forAll (lowerAddr, coeffI)
        {
            // Sign of strong positive upper is opposite of the sign of
            // its diagonal coefficient
            magLower[coeffI] =
                Foam::max(-1*sign(diag[upperAddr[coeffI]])*lower[coeffI], 0);
        }
    }
    else
    {
        magLower = magUpper;
    }

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

        scalarField magScaledDiag = diagFactor_()*magDiag;

        boolList zeroCluster(magDiag.size(), true);

        if (matrix().symmetric())
        {
            // Owner: upper triangle
            forAll (lowerAddr, coeffI)
            {
                if (magUpper[coeffI] > magScaledDiag[lowerAddr[coeffI]])
                {
                    zeroCluster[lowerAddr[coeffI]] = false;
                }
            }

            // Neighbour: lower triangle with symm coefficients
            forAll (upperAddr, coeffI)
            {
                if (magUpper[coeffI] > magScaledDiag[upperAddr[coeffI]])
                {
                    zeroCluster[upperAddr[coeffI]] = false;
                }
            }
        }
        else if (matrix().asymmetric())
        {
            // Owner: upper triangle
            forAll (lowerAddr, coeffI)
            {
                if (magUpper[coeffI] > magScaledDiag[lowerAddr[coeffI]])
                {
                    zeroCluster[lowerAddr[coeffI]] = false;
                }
            }

            // Neighbour: lower triangle with lower coeffs
            forAll (upperAddr, coeffI)
            {
                if (magLower[coeffI] > magScaledDiag[upperAddr[coeffI]])
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

                child_[eqnI] = nCoarseEqns_;
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
        if (child_[rowI] == -1)
        {
            // Found new ungrouped equation
            curEqn = rowI;

            // Reset grouped and upgrouped index
            indexUngrouped = -1;
            indexGrouped = -1;

            // Make next group (coarse equation) and search for neighbours
            child_[curEqn] = nCoarseEqns_;

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

                    curWeight = magUpper[rowCoeffI]/magDiag[curEqn];

                    if (child_[colI] == -1)
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
                    else if (child_[curEqn] != child_[colI])
                    {
                        // Check for neighbour in solo group
                        if (nSolo_ == 0 || child_[colI] != 0)
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

                    curWeight = magLower[losortAddr[rowCoeffI]]/
                        magDiag[curEqn];

                    if (child_[colI] == -1)
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
                    else if (child_[curEqn] != child_[colI])
                    {
                        // Check for neighbour in solo group
                        if (nSolo_ == 0 || child_[colI] != 0)
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

                    child_[nextUngrouped] = child_[curEqn];
                    sizeOfGroups[child_[curEqn]]++;

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
             || sizeOfGroups[child_[nextGrouped]] > maxGroupSize_
            )
            {
                // There is no group to put this equation into
                sizeOfGroups[child_[rowI]]++;
                nCoarseEqns_++;
            }
            else
            {
                // Dump current cell into the best group available
                child_[rowI] = child_[nextGrouped];
                sizeOfGroups[child_[nextGrouped]]++;
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
    if (nCoarseEqns_ > minCoarseEqns() && 3*nCoarseEqns_ <= 2*nRows)
    {
        coarsen_ = true;
    }

    reduce(coarsen_, andOp<bool>());

    if (lduMatrix::debug >= 3)
    {
        // Count solo cells
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clusterAmgPolicy::clusterAmgPolicy
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
    minGroupSize_(groupSize),
    maxGroupSize_(2*groupSize),
    child_(matrix.lduAddr().size()),
    nSolo_(0),
    nCoarseEqns_(0),
    coarsen_(false)
{
    if (groupSize < 2)
    {
        FatalErrorIn
        (
            "clusterAmgPolicy::clusterAmgPolicy\n"
            "(\n"
            "    const lduMatrix& matrix,\n"
            "    const FieldField<Field, scalar>& bouCoeffs,\n"
            "    const FieldField<Field, scalar>& intCoeffs,\n"
            "    const lduInterfaceFieldPtrsList& interfaceFields,\n"
            "    const label groupSize,\n"
            "    const label minCoarseEqns\n"
            ")"
        )   << "Group size smaller than 2 is not allowed"
            << abort(FatalError);
    }

    calcChild();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clusterAmgPolicy::~clusterAmgPolicy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::amgMatrix> Foam::clusterAmgPolicy::restrictMatrix() const
{
    if (!coarsen_)
    {
        FatalErrorIn
        (
            "autoPtr<amgMatrix> clusterAmgPolicy::restrictMatrix() const"
        )   << "Requesting coarse matrix when it cannot be created"
            << abort(FatalError);
    }

    // Construct the coarse matrix and ldu addressing for the next level
    // Algorithm:
    // 1) Loop through all fine coeffs. If the child labels on two sides are
    //    different, this creates a coarse coeff. Define owner and neighbour
    //    for this coeff based on cluster IDs.
    // 2) Check if the coeff has been seen before.  If yes, add the coefficient
    //    to the appropriate field (stored with the equation).  If no, create
    //    a new coeff with neighbour ID and add the coefficient
    // 3) Once all the coeffs have been created, loop through all clusters and
    //    insert the coeffs in the upper order.  At the same time, collect the
    //    owner and neighbour addressing.
    // 4) Agglomerate the diagonal by summing up the fine diagonal

    // Get addressing
    const unallocLabelList& upperAddr = matrix().lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = matrix().lduAddr().lowerAddr();

    const label nFineCoeffs = upperAddr.size();

#   ifdef FULLDEBUG
    if (child_.size() != matrix().lduAddr().size())
    {
        FatalErrorIn
        (
            "autoPtr<amgMatrix> clusterAmgPolicy::restrictMatrix() const"
        )   << "Child array does not correspond to fine level. " << endl
            << " Child size: " << child_.size()
            << " number of equations: " << matrix().lduAddr().size()
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
    // Note: value of coeffRestrictAddr for off-diagonal coefficients
    // touching solo cells will be invalid
    // HJ, 7/Apr/2015
    labelList coeffRestrictAddr(nFineCoeffs);

    // Initial neighbour array (not in upper-triangle order)
    labelList initCoarseNeighb(nFineCoeffs);

    // Counter for coarse coeffs
    label nCoarseCoeffs = 0;

    // Note on zero cluster coarsening
    // If the matrix contains a solo group, it will be in index zero.
    // Since solo equations are disconnected from the rest of the matrix
    // they do not create new off-diagonal coefficients

    // Loop through all fine coeffs
    forAll (upperAddr, fineCoeffI)
    {
        label rmUpperAddr = child_[upperAddr[fineCoeffI]];
        label rmLowerAddr = child_[lowerAddr[fineCoeffI]];

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
            coeffRestrictAddr[fineCoeffI] = -(rmUpperAddr + 1);
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
                    coeffRestrictAddr[fineCoeffI] =
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
                coeffRestrictAddr[fineCoeffI] = nCoarseCoeffs;
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

    forAll (coeffRestrictAddr, fineCoeffI)
    {
        label rmUpperAddr = child_[upperAddr[fineCoeffI]];
        label rmLowerAddr = child_[lowerAddr[fineCoeffI]];

        // If the coefficient touches block zero and solo equations are
        // present, skip it
        if (soloEqns && (rmUpperAddr == 0 || rmLowerAddr == 0))
        {
            continue;
        }

        if (coeffRestrictAddr[fineCoeffI] >= 0)
        {
            coeffRestrictAddr[fineCoeffI] =
                coarseCoeffMap[coeffRestrictAddr[fineCoeffI]];
        }
    }

    // Clear the temporary storage for the coarse matrix data
    blockNnbrs.setSize(0);
    blockNbrsData.setSize(0);
    initCoarseNeighb.setSize(0);
    coarseCoeffMap.setSize(0);


    // Create coarse-level coupled interfaces

    // Create coarse interfaces, addressing and coefficients

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

    // Add the coarse level

    // Set the coarse ldu addressing onto the list
    lduPrimitiveMesh* coarseAddrPtr =
        new lduPrimitiveMesh
        (
            nCoarseEqns_,
            coarseOwner,
            coarseNeighbour,
            true
        );

    // Initialise transfer of restrict addressing on the interface
    forAll (interfaceFields(), intI)
    {
        if (interfaceFields().set(intI))
        {
            const lduInterface& fineInterface =
                interfaceFields()[intI].coupledInterface();

            fineInterface.initInternalFieldTransfer
            (
                Pstream::blocking,
                child_
            );
        }
    }

    // Store neighbour child arrays to avoid tangled communications
    // HJ, 1/Apr/2009
    FieldField<Field, label> fineInterfaceAddr(interfaceFields().size());

    forAll (interfaceFields(), intI)
    {
        if (interfaceFields().set(intI))
        {
            const lduInterface& fineInterface =
                interfaceFields()[intI].coupledInterface();

            fineInterfaceAddr.set
            (
                intI,
                new labelField
                (
                    fineInterface.internalFieldTransfer
                    (
                        Pstream::blocking,
                        child_
                    )
                )
            );
        }
    }

    // Create AMG interfaces
    forAll (interfaceFields(), intI)
    {
        if (interfaceFields().set(intI))
        {
            const lduInterface& fineInterface =
                interfaceFields()[intI].coupledInterface();

            coarseInterfaces.set
            (
                intI,
                AMGInterface::New
                (
                    *coarseAddrPtr,
                    coarseInterfaces,
                    fineInterface,
                    fineInterface.interfaceInternalField(child_),
                    fineInterfaceAddr[intI]
                ).ptr()
            );
        }
    }

    forAll (interfaceFields(), intI)
    {
        if (interfaceFields().set(intI))
        {
            const AMGInterface& coarseInterface =
                refCast<const AMGInterface>(coarseInterfaces[intI]);

            coarseInterfaceFields.set
            (
                intI,
                AMGInterfaceField::New
                (
                    coarseInterface,
                    interfaceFields()[intI]
                ).ptr()
            );

            // Note: scalar agglomeration is done by the interface
            // (always scalar) but in the block matrix it is done by a
            // templated block interface field
            // HJ, 16/Mar/2016
            coarseBouCoeffs.set
            (
                intI,
                coarseInterface.agglomerateCoeffs(bouCoeffs()[intI])
            );

            coarseIntCoeffs.set
            (
                intI,
                coarseInterface.agglomerateCoeffs(intCoeffs()[intI])
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

    // Matrix restriction done!

    // Set the coarse level matrix
    lduMatrix* coarseMatrixPtr = new lduMatrix(*coarseAddrPtr);
    lduMatrix& coarseMatrix = *coarseMatrixPtr;

    // Coarse matrix diagonal initialised by restricting the
    // finer mesh diagonal
    scalarField& coarseDiag = coarseMatrix.diag();
    restrictResidual(matrix().diag(), coarseDiag);

    // Check if matrix is assymetric and if so agglomerate both upper and lower
    // coefficients ...
    if (matrix().hasLower())
    {
        // Get off-diagonal matrix coefficients
        const scalarField& fineUpper = matrix().upper();
        const scalarField& fineLower = matrix().lower();

        // Coarse matrix upper coefficients
        scalarField& coarseUpper = coarseMatrix.upper();
        scalarField& coarseLower = coarseMatrix.lower();

        forAll (coeffRestrictAddr, fineCoeffI)
        {
            label rmUpperAddr = child_[upperAddr[fineCoeffI]];
            label rmLowerAddr = child_[lowerAddr[fineCoeffI]];

            // If the coefficient touches block zero and solo equations are
            // present, skip it
            if (soloEqns && (rmUpperAddr == 0 || rmLowerAddr == 0))
            {
                continue;
            }

            label cCoeff = coeffRestrictAddr[fineCoeffI];

            if (cCoeff >= 0)
            {
                coarseUpper[cCoeff] += fineUpper[fineCoeffI];
                coarseLower[cCoeff] += fineLower[fineCoeffI];
            }
            else
            {
                // Add the fine face coefficients into the diagonal.
                coarseDiag[-1 - cCoeff] +=
                    fineUpper[fineCoeffI] + fineLower[fineCoeffI];
            }
        }
    }
    else // ... Otherwise it is symmetric so agglomerate just the upper
    {
        // Get off-diagonal matrix coefficients
        const scalarField& fineUpper = matrix().upper();

        // Coarse matrix upper coefficients
        scalarField& coarseUpper = coarseMatrix.upper();

        forAll (coeffRestrictAddr, fineCoeffI)
        {
            label rmUpperAddr = child_[upperAddr[fineCoeffI]];
            label rmLowerAddr = child_[lowerAddr[fineCoeffI]];

            // If the coefficient touches block zero and solo equations are
            // present, skip it
            if (soloEqns && (rmUpperAddr == 0 || rmLowerAddr == 0))
            {
                continue;
            }

            label cCoeff = coeffRestrictAddr[fineCoeffI];

            if (cCoeff >= 0)
            {
                coarseUpper[cCoeff] += fineUpper[fineCoeffI];
            }
            else
            {
                // Add the fine face coefficient into the diagonal.
                coarseDiag[-1 - cCoeff] += 2*fineUpper[fineCoeffI];
            }
        }
    }

    // Create and return amgMatrix
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


void Foam::clusterAmgPolicy::restrictResidual
(
    const scalarField& res,
    scalarField& coarseRes
) const
{
    coarseRes = 0;

    forAll (res, i)
    {
        coarseRes[child_[i]] += res[i];
    }
}


void Foam::clusterAmgPolicy::prolongateCorrection
(
    scalarField& x,
    const scalarField& coarseX
) const
{
    forAll (x, i)
    {
        x[i] += coarseX[child_[i]];
    }
}


// ************************************************************************* //
