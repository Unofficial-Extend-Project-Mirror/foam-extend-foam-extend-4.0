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

\*---------------------------------------------------------------------------*/

#include "extendedBlockLduMatrix.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::extendedBlockLduMatrix<Type>::clearOut()
{
    deleteDemandDrivenData(extendedLowerPtr_);
    deleteDemandDrivenData(extendedUpperPtr_);
}


template<class Type>
void Foam::extendedBlockLduMatrix<Type>::mapOffDiagCoeffs
(
    const BlockLduMatrix<Type>& blockLdum
)
{
    if (blockLdum.diagonal())
    {
        WarningIn
        (
            "void extendedBlockLduMatrix<Type>::mapOffDiagCoeffs\n"
            "(\n"
            "    const BlockLduMatrix<Type>& blockLdum\n"
            ")"
        )   << "Attempted to create extended lower/upper coeffs for block "
            << "matrix that is diagonal."
            << nl << endl;
    }
    else if (blockLdum.symmetric())
    {
        // Get reference to faceMap in extended addressing
        const unallocLabelList& faceMap = extLduAddr_.faceMap();

        // Matrix is considered symmetric if the upper is allocated and lower
        // is not allocated. Allocating extended upper only.
        extendedUpperPtr_ = new TypeCoeffField
        (
            extLduAddr_.extendedUpperAddr().size()
        );
        TypeCoeffField& extUpper = *extendedUpperPtr_;

        // Get upper coeffs from underlying lduMatrix
        const TypeCoeffField& upper = blockLdum.upper();

        if (upper.activeType() == blockCoeffBase::SCALAR)
        {
            // Helper type definition
            typedef typename CoeffField<Type>::scalarTypeField activeType;

            // Get references to fields
            const activeType& activeUpper = upper.asScalar();
            activeType& activeExtUpper = extUpper.asScalar();

            // Copy non-zero coeffs from basic lduMatrix into corresponding
            // positions
            forAll (upper, faceI)
            {
                activeExtUpper[faceMap[faceI]] = activeUpper[faceI];
            }
        }
        else if (upper.activeType() == blockCoeffBase::LINEAR)
        {
            // Helper type definition
            typedef typename CoeffField<Type>::linearTypeField activeType;

            // Get references to fields
            const activeType& activeUpper = upper.asLinear();
            activeType& activeExtUpper = extUpper.asLinear();

            // Copy non-zero coeffs from basic lduMatrix into corresponding
            // positions
            forAll (upper, faceI)
            {
                activeExtUpper[faceMap[faceI]] = activeUpper[faceI];
            }
        }
        else if (upper.activeType() == blockCoeffBase::SQUARE)
        {
            // Helper type definition
            typedef typename CoeffField<Type>::squareTypeField activeType;

            // Get references to fields
            const activeType& activeUpper = upper.asSquare();
            activeType& activeExtUpper = extUpper.asSquare();

            // Copy non-zero coeffs from basic lduMatrix into corresponding
            // positions
            forAll (upper, faceI)
            {
                activeExtUpper[faceMap[faceI]] = activeUpper[faceI];
            }
        }
        else
        {
            FatalErrorIn
            (
                "void extendedBlockLduMatrix<Type>::mapOffDiagCoeffs\n"
                "(\n"
                "    const BlockLduMatrix<Type>& blockLdum\n"
                ")"
            )   << "Problem between ordinary block matrix and extended"
                << " block matrix upper coeffs type morphing."
                << abort(FatalError);
        }
    }
    else
    {
        // Get reference to faceMap in extended addressing
        const unallocLabelList& faceMap = extLduAddr_.faceMap();

        // Get number of extended faces
        const label nExtFaces = extLduAddr_.extendedUpperAddr().size();

        // Allocate extended upper and lower
        extendedUpperPtr_ = new TypeCoeffField(nExtFaces);
        TypeCoeffField& extUpper = *extendedUpperPtr_;

        extendedLowerPtr_ = new TypeCoeffField(nExtFaces);
        TypeCoeffField& extLower = *extendedLowerPtr_;

        // Get upper and lower coeffs from underlying lduMatrix
        const TypeCoeffField& upper = blockLdum.upper();
        const TypeCoeffField& lower = blockLdum.lower();

        // Assuming lower and upper have the same active type
        if (upper.activeType() == blockCoeffBase::SCALAR)
        {
            // Helper type definition
            typedef typename CoeffField<Type>::scalarTypeField activeType;

            // Get references to fields
            const activeType& activeUpper = upper.asScalar();
            activeType& activeExtUpper = extUpper.asScalar();
            const activeType& activeLower = lower.asScalar();
            activeType& activeExtLower = extLower.asScalar();

            // Copy non-zero coeffs from basic lduMatrix into corresponding
            // positions
            forAll (upper, faceI)
            {
                activeExtUpper[faceMap[faceI]] = activeUpper[faceI];
                activeExtLower[faceMap[faceI]] = activeLower[faceI];
            }
        }
        else if (upper.activeType() == blockCoeffBase::LINEAR)
        {
            // Helper type definition
            typedef typename CoeffField<Type>::linearTypeField activeType;

            // Get references to fields
            const activeType& activeUpper = upper.asLinear();
            activeType& activeExtUpper = extUpper.asLinear();
            const activeType& activeLower = lower.asLinear();
            activeType& activeExtLower = extLower.asLinear();

            // Copy non-zero coeffs from basic lduMatrix into corresponding
            // positions
            forAll (upper, faceI)
            {
                activeExtUpper[faceMap[faceI]] = activeUpper[faceI];
                activeExtLower[faceMap[faceI]] = activeLower[faceI];
            }
        }
        else if (upper.activeType() == blockCoeffBase::SQUARE)
        {
            // Helper type definition
            typedef typename CoeffField<Type>::squareTypeField activeType;

            // Get references to fields
            const activeType& activeUpper = upper.asSquare();
            activeType& activeExtUpper = extUpper.asSquare();
            const activeType& activeLower = lower.asSquare();
            activeType& activeExtLower = extLower.asSquare();

            // Copy non-zero coeffs from basic lduMatrix into corresponding
            // positions
            forAll (upper, faceI)
            {
                activeExtUpper[faceMap[faceI]] = activeUpper[faceI];
                activeExtLower[faceMap[faceI]] = activeLower[faceI];
            }
        }
        else
        {
            FatalErrorIn
            (
                "void extendedBlockLduMatrix<Type>::mapOffDiagCoeffs\n"
                "(\n"
                "    const BlockLduMatrix<Type>& blockLdum\n"
                ")"
            )   << "Problem between ordinary block matrix and extended"
                << " block matrix upper/lower coeffs type morphing."
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::extendedBlockLduMatrix<Type>::extendedBlockLduMatrix
(
    const BlockLduMatrix<Type>& blockLdum,
    const extendedLduAddressing& extLduAddr
)
:
    basicBlockLduMatrix_(blockLdum),
    extLduAddr_(extLduAddr),
    extendedLowerPtr_(NULL),
    extendedUpperPtr_(NULL)
{
    if (debug)
    {
        InfoIn
        (
            "extendedBlockLduMatrix<Type>::extendedBlockLduMatrix\n"
            "(\n"
            "    const BlockLduMatrix<Type>& blockLdum,\n"
            "    const extendedLduAddressing& extLduAddr\n"
            ")"
        )   << "Constructing extendedBlockLduMatrix."
            << endl;
    }

    // Map off diag coeffs from original block matrix to this extended block
    // matrix
    mapOffDiagCoeffs(blockLdum);
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

template<class Type>
Foam::extendedBlockLduMatrix<Type>::~extendedBlockLduMatrix()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
typename Foam::extendedBlockLduMatrix<Type>::TypeCoeffField&
Foam::extendedBlockLduMatrix<Type>::extendedLower()
{
    if (!extendedLowerPtr_)
    {
        if (extendedUpperPtr_)
        {
            extendedLowerPtr_ = new TypeCoeffField(*extendedUpperPtr_);
        }
        else
        {
            extendedLowerPtr_ = new TypeCoeffField
            (
                extLduAddr_.extendedLowerAddr().size()
            );
        }
    }

    return *extendedLowerPtr_;
}


template<class Type>
typename Foam::extendedBlockLduMatrix<Type>::TypeCoeffField&
Foam::extendedBlockLduMatrix<Type>::extendedUpper()
{
    if (!extendedUpperPtr_)
    {
        if (extendedLowerPtr_)
        {
            extendedUpperPtr_ = new TypeCoeffField(*extendedLowerPtr_);
        }
        else
        {
            extendedUpperPtr_ = new TypeCoeffField
            (
                extLduAddr_.extendedUpperAddr().size()
            );
        }
    }

    return *extendedUpperPtr_;
}


template<class Type>
const typename Foam::extendedBlockLduMatrix<Type>::TypeCoeffField&
Foam::extendedBlockLduMatrix<Type>::extendedLower() const
{
    if (!extendedLowerPtr_ && !extendedUpperPtr_)
    {
        FatalErrorIn
        (
            "const CoeffField<Type>& "
            "extendedBlockLduMatrix<Type>::extendedLower() const"
        )   << "extendedLowerPtr_ or extendedUpperPtr_ unallocated"
            << abort(FatalError);
    }

    if (extendedLowerPtr_)
    {
        return *extendedLowerPtr_;
    }
    else
    {
        return *extendedUpperPtr_;
    }
}


template<class Type>
const typename Foam::extendedBlockLduMatrix<Type>::TypeCoeffField&
Foam::extendedBlockLduMatrix<Type>::extendedUpper() const
{
    if (!extendedLowerPtr_ && !extendedUpperPtr_)
    {
        FatalErrorIn
        (
            "const CoeffField<Type>& "
            "extendedBlockLduMatrix<Type>::extendedLower() const"
        )   << "extendedLowerPtr_ or extendedUpperPtr_ unallocated"
            << abort(FatalError);
    }

    if (extendedUpperPtr_)
    {
        return *extendedUpperPtr_;
    }
    else
    {
        return *extendedLowerPtr_;
    }
}


// ************************************************************************* //
