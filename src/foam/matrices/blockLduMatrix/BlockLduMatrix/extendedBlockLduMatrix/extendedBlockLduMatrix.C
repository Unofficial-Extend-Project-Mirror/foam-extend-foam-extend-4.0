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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::extendedBlockLduMatrix<Type>::extendedBlockLduMatrix
(
    const blockLduMatrix<Type>& blockLdum,
    const label extensionLevel,
    const polyMesh& polyMesh
)
:
    basicBlockLduMatrix_(blockLdum),
    extLduAddr_
    (
        extendedLduAddressing::New
        (
            polyMesh,
            blockLdum.lduAddr(),
            extensionLevel
        )
    ),
    extendedLowerPtr_(NULL),
    extendedUpperPtr_(NULL)
{
    if (debug)
    {
        Info<< "extendedBlockLduMatrix(lduMatrix&, label, polyMesh&) :"
               "Constructing extendedBlockLduMatrix."
            << endl;
    }

    if (blockLdum.diagonal())
    {
        WarningIn("extendedBlockLduMatrix(lduMatrix&, label, polyMesh&)")
            << "Attempted to create extended lower/upper coeffs for block "
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

            // Get upper coeffs from underlying lduMatrix
            const TypeCoeffField& upper = blockLdum.upper();

            // Copy non-zero coeffs from basic lduMatrix into corresponding
            // positions
            forAll (upper, faceI)
            {
                extUpper[faceMap[faceI]] = upper[faceI];
            }
        }
        else
        {
            // Allocate extended lower only
            extendedLowerPtr_ = new TypeCoeffField
            (
                extLduAddr_.extendedLowerAddr().size(),
                pTraits<Type>::zero
            );
            TypeCoeffField& extLower = *extendedLowerPtr_;

            // Get lower coeffs from underlying lduMatrix
            const TypeCoeffField& lower = blockLdum.lower();

            // Copy non-zero coeffs from basic lduMatrix into corresponding
            // positions
            forAll (lower, faceI)
            {
                extLower[faceMap[faceI]] = lower[faceI];
            }
        }
    }
    else
    {
        // Get reference to faceMap in extended addressing
        const unallocLabelList& faceMap = extLduAddr_.faceMap();

        // Get number of extended faces
        const label nExtFaces = extLduAddr_.extendedUpperAddr().size();

        // Allocate extended upper and lower
        extendedUpperPtr_ = new TypeCoeffField(nExtFaces, pTraits<Type>::zero);
        TypeCoeffField& extUpper = *extendedUpperPtr_;

        extendedLowerPtr_ = new TypeCoeffField(nExtFaces, pTraits<Type>::zero);
        TypeCoeffField& extLower = *extendedLowerPtr_;

        // Get upper and lower coeffs from underlying lduMatrix
        const TypeCoeffField& upper = blockLdum.upper();
        const TypeCoeffField& lower = blockLdum.lower();

        // Copy non-zero coeffs from basic lduMatrix into corresponding
        // positions
        forAll (upper, faceI)
        {
            extUpper[faceMap[faceI]] = upper[faceI];
            extLower[faceMap[faceI]] = lower[faceI];
        }
    }
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
                extLduAddr_.extendedLowerAddr().size(),
                pTraits<Type>::zero
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
                extLduAddr_.extendedUpperAddr().size(),
                pTraits<Type>::zero
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
