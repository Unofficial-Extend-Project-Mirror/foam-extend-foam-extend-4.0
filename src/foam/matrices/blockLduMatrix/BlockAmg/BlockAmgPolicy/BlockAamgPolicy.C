/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    BlockAamgPolicy

Description
    Agglomerative AMG policy, adjusted for BlockLduMatrix

Author
    Klas Jareteg, 2012-12-13

\*---------------------------------------------------------------------------*/

#include "BlockAamgPolicy.H"
#include "coeffFields.H"
#include "addToRunTimeSelectionTable.H"
#include "BlockGAMGInterfaceField.H"
#include "processorLduInterfaceField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::scalar Foam::BlockAamgPolicy<Type>::weightFactor_ = 0.65;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::BlockAamgPolicy<Type>::calcChild()
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
        // Use child array to count number of entries per row.
        // Reset the list for counting
        labelList& nPerRow = child_;
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

        // Reset child array
        child_ = -1;
    }


    // Calculate agglomeration

    // Get matrix coefficients
    const CoeffField<Type>& diag = matrix_.diag();
    const CoeffField<Type>& upper = matrix_.upper();

    labelList sizeOfGroups(nRows, 0);

    nCoarseEqns_ = 0;

    label indexU, indexG, colI, curEqn, nextEqn, groupPassI;

    scalar magRowDiag, magColDiag, weight, weightU, weightG;

    // Magnitudes pre-calculated
    Field<scalar> magDiag(diag.size());
    Field<scalar> magUpper(upper.size());

    normPtr_->coeffMag(diag,magDiag);
    normPtr_->coeffMag(upper,magUpper);

    for (label eqnI = 0; eqnI < nRows; eqnI++)
    {
        if (child_[eqnI] == -1)
        {
            curEqn = eqnI;

            indexU = -1;
            indexG = -1;

            child_[curEqn] = nCoarseEqns_;

            magRowDiag = magDiag[curEqn];

            for (groupPassI = 1; groupPassI < groupSize_; groupPassI++)
            {
                weightU = 0;
                weightG = 0;

                indexU = -1;
                indexG = -1;

                for
                (
                    label rowCoeffI = rowOffset[curEqn];
                    rowCoeffI < rowOffset[curEqn + 1];
                    rowCoeffI++
                )
                {
                    colI = cols[rowCoeffI];

                    magColDiag = magDiag[colI];

                    weight = (magUpper[cIndex[rowCoeffI]]
                        / max(magRowDiag, magColDiag));

                    if (child_[colI] == -1)
                    {
                        if (indexU == -1 || weight > weightU)
                        {
                            indexU = rowCoeffI;
                            weightU = weight;
                        }
                    }
                    else if (child_[curEqn] != child_[colI])
                    {
                        if (indexG == -1 || weight > weightG)
                        {
                            indexG = rowCoeffI;
                            weightG = weight;
                        }
                    }
                }

                if
                (
                    indexU != -1
                 && (indexG == -1 || weightU >= weightFactor_*weightG)
                )
                {
                    // Found new element of group.  Add it and use as
                    // start of next search

                    nextEqn = cols[indexU];

                    child_[nextEqn] = child_[curEqn];
                    sizeOfGroups[child_[curEqn]];

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
             || indexG == -1
             || sizeOfGroups[child_[cols[indexG]]] > (groupSize_ + 2)
            )
            {
                sizeOfGroups[child_[eqnI]]++;
                nCoarseEqns_++;
            }
            else
            {
                child_[eqnI] = child_[cols[indexG]];
                sizeOfGroups[child_[cols[indexG]]]++;
            }
        }
    }

    // The decision on parallel agglomeration needs to be made for the
    // whole gang of processes; otherwise I may end up with a different
    // number of agglomeration levels on different processors.

    if
    (
        nCoarseEqns_ > BlockAmgPolicy<Type>::minCoarseEqns()
     && 3*nCoarseEqns_ <= 2*nRows
    )
    {
        coarsen_ = true;
    }

    reduce(coarsen_, andOp<bool>());

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockAamgPolicy<Type>::BlockAamgPolicy
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict,
    const label groupSize,
    const label minCoarseEqns
)
:
    BlockAmgPolicy<Type>(matrix, dict, groupSize, minCoarseEqns),
    matrix_(matrix),
    normPtr_
    (
        BlockCoeffNorm<Type>::New
        (
            dict
        )
    ),
    child_(matrix_.lduAddr().size()),
    groupSize_(groupSize),
    nCoarseEqns_(0),
    coarsen_(false)
{
    calcChild();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockAamgPolicy<Type>::~BlockAamgPolicy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockLduMatrix<Type> >
Foam::BlockAamgPolicy<Type>::restrictMatrix() const
{
    if (!coarsen_)
    {
        FatalErrorIn
        (
            "autoPtr<amgMatrix> BlockAamgPolicy<Type>::restrictMatrix() const"
        )   << "Requesting coarse matrix when it cannot be created"
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

    const label nFineCoeffs = upperAddr.size();

#   ifdef FULLDEBUG
    if (child_.size() != matrix_.lduAddr().size())
    {
        FatalErrorIn
        (
            "autoPtr<BlockLduMatrix<Type> > BlockAamgPolicy<Type>::"
            "restrictMatrix() const"
        )   << "Child array does not correspond to fine level. " << endl
            << " Child size: " << child_.size()
            << " number of equations: " << matrix_.lduAddr().size()
            << abort(FatalError);
    }
#   endif


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
        label rmUpperAddr = child_[upperAddr[fineCoeffi]];
        label rmLowerAddr = child_[lowerAddr[fineCoeffi]];

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
    lduInterfacePtrsList* coarseInterfacesPtr =
        new lduInterfacePtrsList(interfaceSize);
    lduInterfacePtrsList& coarseInterfaces = *coarseInterfacesPtr;

    labelListList coarseInterfaceAddr(interfaceSize);

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
            interfaceFields[intI].interface().initInternalFieldTransfer
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
                interfaceFields[intI].interface();

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

    // Create GAMG interfaces
    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            const lduInterface& fineInterface =
                interfaceFields[intI].interface();

            coarseInterfaces.set
            (
                intI,
                GAMGInterface::New
                (
                    *coarseAddrPtr,
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
            const GAMGInterface& coarseInterface =
                refCast<const GAMGInterface>(coarseInterfaces[intI]);

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

    // Set the coarse level matrix
    BlockLduMatrix<Type>* coarseMatrixPtr =
        new BlockLduMatrix<Type>(*coarseAddrPtr);
    BlockLduMatrix<Type>& coarseMatrix = *coarseMatrixPtr;

    typename BlockLduInterfaceFieldPtrsList<Type>::Type&
        coarseInterfaceFieldsTransfer =
          coarseMatrix.interfaces();

    // Aggolmerate the upper and lower coupled coefficients
    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            const GAMGInterface& coarseInterface =
                refCast<const GAMGInterface>(coarseInterfaces[intI]);

            coarseInterfaceFieldsTransfer.set
            (
                intI,
                BlockGAMGInterfaceField<Type>::New
                (
                    coarseInterface,
                    (interfaceFields[intI])
                ).ptr()
            );

            coarseMatrix.coupleUpper().set
            (
                intI,
                coarseInterface.agglomerateBlockCoeffs
                (
                    matrix_.coupleUpper()[intI]
                )
            );

            coarseMatrix.coupleLower().set
            (
                intI,
                coarseInterface.agglomerateBlockCoeffs
                (
                    matrix_.coupleLower()[intI]
                )
            );
        }
    }

    // Matrix restriction done!

    typedef CoeffField<Type> TypeCoeffField;

    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    TypeCoeffField& coarseUpper = coarseMatrix.upper();
    TypeCoeffField& coarseDiag = coarseMatrix.diag();
    const TypeCoeffField& fineUpper = matrix_.upper();
    const TypeCoeffField& fineDiag = matrix_.diag();

    // KRJ: 2013-01-31: Many cases needed as there are different combinations
    if (matrix_.asymmetric())
    {
        TypeCoeffField& coarseLower = coarseMatrix.lower();
        const TypeCoeffField& fineLower = matrix_.lower();

        if (fineDiag.activeType() == blockCoeffBase::SQUARE)
        {
            if (fineUpper.activeType() == blockCoeffBase::SQUARE)
            {
                squareTypeField& activeCoarseUpper = coarseUpper.asSquare();
                squareTypeField& activeCoarseDiag = coarseDiag.asSquare();
                const squareTypeField& activeFineUpper = fineUpper.asSquare();

                squareTypeField& activeCoarseLower = coarseLower.asSquare();
                const squareTypeField& activeFineLower = fineLower.asSquare();

                restrictResidual(fineDiag, coarseDiag);

                forAll(coeffRestrictAddr, fineCoeffI)
                {
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
                            activeFineUpper[fineCoeffI] + activeFineLower[fineCoeffI];
                    }
                }
            }
            else if (fineUpper.activeType() == blockCoeffBase::LINEAR)
            {
                FatalErrorIn("autoPtr<amgMatrix> BlockAamgPolicy<Type>::restrictMatrix() const")
                    << "Matrix diagonal of square type and upper of linear type is not implemented"
                    << abort(FatalError);
            }
            else
            {
                FatalErrorIn("autoPtr<amgMatrix> BlockAamgPolicy<Type>::restrictMatrix() const")
                    << "Matrix diagonal of square type and upper of scalar type is not implemented"
                    << abort(FatalError);
            }
        }
        else if (fineDiag.activeType() == blockCoeffBase::LINEAR)
        {
                FatalErrorIn("autoPtr<amgMatrix> BlockAamgPolicy<Type>::restrictMatrix() const")
                    << "Matrix diagonal of linear type not implemented"
                    << abort(FatalError);
        }
        else
        {
                FatalErrorIn("autoPtr<amgMatrix> BlockAamgPolicy<Type>::restrictMatrix() const")
                    << "Matrix diagonal of scalar type not implemented"
                    << abort(FatalError);
        }
    }
    else
    {
        if (fineDiag.activeType() == blockCoeffBase::SQUARE)
        {
            if (fineUpper.activeType() == blockCoeffBase::SQUARE)
            {
                squareTypeField& activeCoarseUpper = coarseUpper.asSquare();
                squareTypeField& activeCoarseDiag = coarseDiag.asSquare();
                const squareTypeField& activeFineUpper = fineUpper.asSquare();

                restrictResidual(fineDiag, coarseDiag);

                forAll(coeffRestrictAddr, fineCoeffI)
                {
                    label cCoeff = coeffRestrictAddr[fineCoeffI];

                    if (cCoeff >= 0)
                    {
                        activeCoarseUpper[cCoeff] += activeFineUpper[fineCoeffI];
                    }
                    else
                    {
                        // Add the fine face coefficient into the diagonal.
                        activeCoarseDiag[-1 - cCoeff] += 2*activeFineUpper[fineCoeffI];
                    }
                }
            }
            else if (fineUpper.activeType() == blockCoeffBase::LINEAR)
            {
                FatalErrorIn("autoPtr<amgMatrix> BlockAamgPolicy<Type>::restrictMatrix() const")
                    << "Matrix diagonal of square type and upper of linear type is not implemented"
                    << abort(FatalError);
            }
            else
            {
                FatalErrorIn("autoPtr<amgMatrix> BlockAamgPolicy<Type>::restrictMatrix() const")
                    << "Matrix diagonal of square type and upper of scalar type is not implemented"
                    << abort(FatalError);
            }

        }
        else if (fineDiag.activeType() == blockCoeffBase::LINEAR)
        {
                FatalErrorIn
                (
                    "autoPtr<amgMatrix> BlockAamgPolicy<Type>::"
                    "restrictMatrix() const"
                )   << "Matrix diagonal of linear type not implemented"
                    << abort(FatalError);
        }
        else
        {
                FatalErrorIn
                (
                    "autoPtr<amgMatrix> BlockAamgPolicy<Type>::"
                    "restrictMatrix() const"
                )   << "Matrix diagonal of scalar type not implemented"
                    << abort(FatalError);
        }
    }

    return autoPtr<BlockLduMatrix<Type> >
    (
        new BlockLduMatrix<Type>
        (
            coarseMatrix
        )
    );
}


template<class Type>
void Foam::BlockAamgPolicy<Type>::restrictResidual
(
    const CoeffField<Type>& res,
    CoeffField<Type>& coarseRes
) const
{
    typedef CoeffField<Type> TypeCoeffField;

    if (res.activeType() == blockCoeffBase::SQUARE &&
            coarseRes.activeType() == blockCoeffBase::SQUARE)
    {
        typedef typename TypeCoeffField::squareTypeField squareTypeField;

        squareTypeField& activeCoarseRes = coarseRes.asSquare();
        const squareTypeField& activeRes = res.asSquare();

        // KRJ: 2013-02-05: To be done by correct pTraits?
        forAll (coarseRes, i)
        {
            activeCoarseRes[i] *= 0.0;
        }

        forAll (res, i)
        {
            activeCoarseRes[child_[i]] += activeRes[i];
        }
    }
    else
    {
        FatalErrorIn("void  BlockAamgPolicy<Type>::restrictResidual() const")
            << "Only present for square type coeff type"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::BlockAamgPolicy<Type>::restrictResidual
(
    const Field<Type>& res,
    Field<Type>& coarseRes
) const
{
    coarseRes = pTraits<Type>::zero;

    forAll (res, i)
    {
        coarseRes[child_[i]] += res[i];
    }
}


template<class Type>
void Foam::BlockAamgPolicy<Type>::prolongateCorrection
(
    Field<Type>& x,
    const Field<Type>& coarseX
) const
{
    forAll (x, i)
    {
        x[i] += coarseX[child_[i]];
    }
}


// ************************************************************************* //
