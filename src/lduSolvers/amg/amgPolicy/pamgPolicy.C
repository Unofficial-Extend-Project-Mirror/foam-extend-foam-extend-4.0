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
    pamgPolicy

Description
    Pairwise agglomerative AMG policy.  Legacy implementation

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "pamgPolicy.H"
#include "amgMatrix.H"
#include "boolList.H"
#include "addToRunTimeSelectionTable.H"
#include "AMGInterfaceField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pamgPolicy, 0);

    addToRunTimeSelectionTable(amgPolicy, pamgPolicy, matrix);

} // End namespace Foam


const Foam::debug::tolerancesSwitch
Foam::pamgPolicy::diagFactor_
(
    "pamgDiagFactor",
    1e-8
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pamgPolicy::calcChild()
{
    // Algorithm:
    // 1) Create temporary equation addressing using a double-pass algorithm.
    //    to create the offset table.
    // 2) Loop through all equations and for each equation find the best fit
    //    neighbour.  If all neighbours are grouped, add equation to best group

    // Get addressing
    const label nEqns = matrix_.lduAddr().size();

    const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();

    // Get off-diagonal matrix coefficients
    const scalarField& upper = matrix_.upper();

    // For each equation calculate coeffs
    labelList cols(upperAddr.size() + lowerAddr.size());
    labelList rowOffsets(nEqns + 1);

    // Memory management
    {
        labelList nNbrs(nEqns, 0);

        forAll (upperAddr, coeffI)
        {
            nNbrs[upperAddr[coeffI]]++;
        }

        forAll (lowerAddr, coeffI)
        {
            nNbrs[lowerAddr[coeffI]]++;
        }

        rowOffsets[0] = 0;

        forAll (nNbrs, eqnI)
        {
            rowOffsets[eqnI + 1] = rowOffsets[eqnI] + nNbrs[eqnI];
        }

        // Reset the list to use as counter
        nNbrs = 0;

        forAll (upperAddr, coeffI)
        {
            cols
            [
                rowOffsets[upperAddr[coeffI]] + nNbrs[upperAddr[coeffI]]
            ] = coeffI;

            nNbrs[upperAddr[coeffI]]++;
        }

        forAll (lowerAddr, coeffI)
        {
            cols
            [
                rowOffsets[lowerAddr[coeffI]] + nNbrs[lowerAddr[coeffI]]
            ] = coeffI;

            nNbrs[lowerAddr[coeffI]]++;
        }
    }

    // Reset child array
    child_.setSize(nEqns, -1);

    // Calculate agglomeration

    // Get matrix coefficients
    const scalarField& diag = matrix_.diag();

    scalarField magOffDiag;

    if (matrix_.asymmetric())
    {
        magOffDiag = Foam::max(mag(matrix_.upper()), mag(matrix_.lower()));
    }
    else if (matrix_.symmetric())
    {
        magOffDiag = mag(matrix_.upper());
    }
    else
    {
        // Diag only matrix.  Reset and return
        child_ = 0;
        nCoarseEqns_ = 1;

        return;
    }

    nCoarseEqns_ = 0;

    // Gather disconnected and weakly connected equations into cluster zero
    // Weak connection is assumed to be the one where the off-diagonal
    // coefficient is smaller than diagFactor_*diag
    {
        // Algorithm
        // Mark all cells to belong to zero cluster
        // Go through all upper and lower coefficients and for the ones
        // larger than threshold mark the equations out of cluster zero

        scalarField magScaledDiag = diagFactor_()*mag(diag);

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

                child_[eqnI] = nCoarseEqns_;
            }
        }

        if (nSolo_ > 0)
        {
            // Found solo equations
            nCoarseEqns_++;

            if (lduMatrix::debug >= 2)
            {
                Pout<< "Found " << nSolo_ << " weakly connected equations."
                    << endl;
            }
        }
    }

    // Go through the off-diagonal and create clusters, marking the child array
    for (label eqnI = 0; eqnI < nEqns; eqnI++)
    {
        if (child_[eqnI] < 0)
        {
            label matchCoeffNo = -1;
            scalar maxCoeff = -GREAT;

            // Check row to find ungrouped neighbour with largest coefficient
            for
            (
                label rowCoeffI = rowOffsets[eqnI];
                rowCoeffI < rowOffsets[eqnI + 1];
                rowCoeffI++
            )
            {
                label coeffI = cols[rowCoeffI];

                // I don't know whether the current equation is owner
                //  or neighbour.  Therefore I'll check both sides
                if
                (
                    child_[upperAddr[coeffI]] < 0
                 && child_[lowerAddr[coeffI]] < 0
                 && mag(upper[coeffI]) > maxCoeff
                )
                {
                    // Match found. Pick up all the necessary data
                    matchCoeffNo = coeffI;
                    maxCoeff = mag(upper[coeffI]);
                }
            }

            if (matchCoeffNo >= 0)
            {
                // Make a new group
                child_[upperAddr[matchCoeffNo]] = nCoarseEqns_;
                child_[lowerAddr[matchCoeffNo]] = nCoarseEqns_;
                nCoarseEqns_++;
            }
            else
            {
                // No match. Find the best neighbouring cluster and
                // put the equation there
                label clusterMatchCoeffNo = -1;
                scalar clusterMaxCoeff = -GREAT;

                for
                (
                    label rowCoeffI = rowOffsets[eqnI];
                    rowCoeffI < rowOffsets[eqnI + 1];
                    rowCoeffI++
                )
                {
                    label coeffI = cols[rowCoeffI];

                    if (mag(upper[coeffI]) > clusterMaxCoeff)
                    {
                        clusterMatchCoeffNo = coeffI;
                        clusterMaxCoeff = mag(upper[coeffI]);
                    }
                }

                if (clusterMatchCoeffNo >= 0)
                {
                    // Add the equation to the best cluster
                    child_[eqnI] = max
                    (
                        child_[upperAddr[clusterMatchCoeffNo]],
                        child_[lowerAddr[clusterMatchCoeffNo]]
                    );
                }
                else
                {
                    // This is a Dirichlet point: no neighbours
                    // Put it into its own cluster
                    child_[eqnI] = nCoarseEqns_;
                    nCoarseEqns_++;
                }
            }
        }
    }

    // Reverse the map to facilitate later agglomeration.
    // Keep zero cluster at zero index
    forAll (child_, eqnI)
    {
        child_[eqnI] = nCoarseEqns_ - 1 - child_[eqnI];
    }

    // The decision on parallel agglomeration needs to be made for the
    // whole gang of processes; otherwise I may end up with a different
    // number of agglomeration levels on different processors.

    if (nCoarseEqns_ > minCoarseEqns() && 3*nCoarseEqns_ <= 2*nEqns)
    {
        coarsen_ = true;
    }

    reduce(coarsen_, andOp<bool>());

    if (lduMatrix::debug >= 2)
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
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pamgPolicy::pamgPolicy
(
    const lduMatrix& matrix,
    const label groupSize,
    const label minCoarseEqns
)
:
    amgPolicy(groupSize, minCoarseEqns),
    matrix_(matrix),
    child_(),
    nSolo_(0),
    nCoarseEqns_(0),
    coarsen_(false)
{
    calcChild();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pamgPolicy::~pamgPolicy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::amgMatrix> Foam::pamgPolicy::restrictMatrix
(
    const FieldField<Field, scalar>& bouCoeffs,
    const FieldField<Field, scalar>& intCoeffs,
    const lduInterfaceFieldPtrsList& interfaceFields
) const
{
    if (!coarsen_)
    {
        FatalErrorIn("autoPtr<amgMatrix> pamgPolicy::restrictMatrix() const")
            << "Requesting coarse matrix when it cannot be created"
            << abort(FatalError);
    }

    // Construct the coarse matrix and ldu addressing for the next level
    // Algorithm:
    // 1) Loop through all fine coeffs. If the child labels on two sides are
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

    label nFineCoeffs = upperAddr.size();

#   ifdef FULLDEBUG
    if (child_.size() != matrix_.lduAddr().size())
    {
        FatalErrorIn
        (
            "autoPtr<amgMatrix> pamgPolicy::restrictMatrix() const"
        )   << "Child array does not correspond to fine level. " << endl
            << " Child size: " << child_.size()
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

    scalarList blockCoeffsData(maxNnbrs*nCoarseEqns_, 0.0);

    // Create face-restriction addressing
    // Note: value of coeffRestrictAddr for off-diagonal coefficients
    // touching solo cells will be invalid
    // HJ, 7/Apr/2015
    labelList coeffRestrictAddr(nFineCoeffs);

    // Initial neighbour array (not in upper-triangle order)
    labelList initCoarseNeighb(nFineCoeffs);

    // Counter for coarse coeffs
    label nCoarseCoeffs = 0;

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
            label* ccCoeffs = &blockNbrsData[maxNnbrs*cOwn];

            bool nbrFound = false;
            label& ccnCoeffs = blockNnbrs[cOwn];

            for (int i = 0; i < ccnCoeffs; i++)
            {
                if (initCoarseNeighb[ccCoeffs[i]] == cNei)
                {
                    nbrFound = true;
                    coeffRestrictAddr[fineCoeffI] = ccCoeffs[i];
                    break;
                }
            }

            if (!nbrFound)
            {
                if (ccnCoeffs >= maxNnbrs)
                {
                    label oldMaxNnbrs = maxNnbrs;
                    maxNnbrs *= 2;

                    blockNbrsData.setSize(maxNnbrs*nCoarseEqns_);

                    forAllReverse(blockNnbrs, i)
                    {
                        label* oldCcNbrs = &blockNbrsData[oldMaxNnbrs*i];
                        label* newCcNbrs = &blockNbrsData[maxNnbrs*i];

                        for (int j=0; j<blockNnbrs[i]; j++)
                        {
                            newCcNbrs[j] = oldCcNbrs[j];
                        }
                    }

                    ccCoeffs = &blockNbrsData[maxNnbrs*cOwn];
                }

                ccCoeffs[ccnCoeffs] = nCoarseCoeffs;
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
        new lduInterfacePtrsList(interfaceFields.size());
    lduInterfacePtrsList& coarseInterfaces = *coarseInterfacesPtr;

    // Set the coarse interfaceFields and coefficients
    lduInterfaceFieldPtrsList* coarseInterfaceFieldsPtr =
        new lduInterfaceFieldPtrsList(interfaceFields.size());
    lduInterfaceFieldPtrsList& coarseInterfaceFields =
        *coarseInterfaceFieldsPtr;

    FieldField<Field, scalar>* coarseBouCoeffsPtr =
        new FieldField<Field, scalar>(interfaceFields.size());
    FieldField<Field, scalar>& coarseBouCoeffs = *coarseBouCoeffsPtr;

    FieldField<Field, scalar>* coarseIntCoeffsPtr =
        new FieldField<Field, scalar>(interfaceFields.size());
    FieldField<Field, scalar>& coarseIntCoeffs = *coarseIntCoeffsPtr;

    labelListList coarseInterfaceAddr(interfaceFields.size());

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
    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            interfaceFields[intI].coupledInterface().initInternalFieldTransfer
            (
                Pstream::blocking,
                child_
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
                        child_
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
                    *coarseAddrPtr,
                    coarseInterfaces,
                    fineInterface,
                    fineInterface.interfaceInternalField(child_),
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

            coarseInterfaceFields.set
            (
                intI,
                AMGInterfaceField::New
                (
                    coarseInterface,
                    interfaceFields[intI]
                ).ptr()
            );

            // Note: scalar agglomeration is done by the interface
            // (always scalar) but in the block matrix it is done by a
            // templated block interface field
            // HJ, 16/Mar/2016
            coarseBouCoeffs.set
            (
                intI,
                coarseInterface.agglomerateCoeffs(bouCoeffs[intI])
            );

            coarseIntCoeffs.set
            (
                intI,
                coarseInterface.agglomerateCoeffs(intCoeffs[intI])
            );

            coarseInterfaceAddr[intI] = coarseInterface.faceCells();
        }
    }

    // Add interfaces
    coarseAddrPtr->addInterfaces
    (
        *coarseInterfacesPtr,
        coarseInterfaceAddr,
        matrix_.patchSchedule()
    );

    // Matrix restriction done!

    // Set the coarse level matrix
    lduMatrix* coarseMatrixPtr = new lduMatrix(*coarseAddrPtr);
    lduMatrix& coarseMatrix = *coarseMatrixPtr;

    // Coarse matrix diagonal initialised by restricting the
    // finer mesh diagonal
    scalarField& coarseDiag = coarseMatrix.diag();
    restrictResidual(matrix_.diag(), coarseDiag);

    // Check if matrix is assymetric and if so agglomerate both upper and lower
    // coefficients ...
    if (matrix_.hasLower())
    {
        // Get off-diagonal matrix coefficients
        const scalarField& fineUpper = matrix_.upper();
        const scalarField& fineLower = matrix_.lower();

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
        const scalarField& fineUpper = matrix_.upper();

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


void Foam::pamgPolicy::restrictResidual
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


void Foam::pamgPolicy::prolongateCorrection
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
