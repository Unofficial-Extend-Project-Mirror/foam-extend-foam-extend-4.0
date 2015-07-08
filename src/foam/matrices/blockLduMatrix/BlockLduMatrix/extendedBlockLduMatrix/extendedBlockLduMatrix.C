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

\*---------------------------------------------------------------------------*/

#include "extendedBlockLduMatrix.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::extendedBlockLduMatrix<Type>::mapOffDiagCoeffs
(
    const blockLduMatrix<Type>& blockLdum,
    const label extensionLevel,
    const polyMesh& polyMesh
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

        // Avoid assuming it's upper if the matrix is symmetric
        if (blockLdum.hasUpper())
        {
            // Allocate extended upper only
            extendedUpperPtr_ = new TypeCoeffField
            (
                extLduAddr_.extendedUpperAddr().size(),
                pTraits<Type>::zero
            );
            TypeCoeffField& extUpper = *extendedUpperPtr_;

            // Get references to fields
            const activeType& activeUpper = upper.asScalar();
            activeType& activeExtUpper = extUpper.asScalar();

            // Copy non-zero coeffs from basic lduMatrix into corresponding
            // positions
            forAll (upper, faceI)
            {
                extUpper[faceMap[faceI]] = upper[faceI];
            }
        }
        else if (upper.activeType() == blockCoeffBase::LINEAR)
        {
            // Allocate extended lower only
            extendedLowerPtr_ = new TypeCoeffField
            (
                extLduAddr_.extendedLowerAddr().size(),
                pTraits<Type>::zero
            );
            TypeCoeffField& extLower = *extendedLowerPtr_;

            // Get references to fields
            const activeType& activeUpper = upper.asLinear();
            activeType& activeExtUpper = extUpper.asLinear();

            // Copy non-zero coeffs from basic lduMatrix into corresponding
            // positions
            forAll (lower, faceI)
            {
                extLower[faceMap[faceI]] = lower[faceI];
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
            extUpper[faceMap[faceI]] = upper[faceI];
            extLower[faceMap[faceI]] = lower[faceI];
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
    if (extendedLowerPtr_)
    {
        delete extendedLowerPtr_;
    }

    if (extendedUpperPtr_)
    {
        delete extendedUpperPtr_;
    }
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
