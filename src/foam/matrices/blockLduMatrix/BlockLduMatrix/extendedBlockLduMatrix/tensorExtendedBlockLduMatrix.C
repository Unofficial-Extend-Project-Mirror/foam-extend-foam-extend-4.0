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

#ifndef tensorExtendedBlockLduMatrix_H
#define tensorExtendedBlockLduMatrix_H

#include "extendedBlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void Foam::extendedBlockLduMatrix<Foam::tensor>::mapOffDiagCoeffs
(
    const BlockLduMatrix<tensor>& blockLdum
)
{
    if (blockLdum.diagonal())
    {
        WarningIn
        (
            "void extendedBlockLduMatrix<tensor>::mapOffDiagCoeffs\n"
            "(\n"
            "    const BlockLduMatrix<tensor>& blockLdum\n"
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
        if (blockLdum.thereIsUpper())
        {
            // Allocate extended upper only
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
                typedef CoeffField<tensor>::scalarTypeField
                    activeType;

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
                typedef CoeffField<tensor>::linearTypeField
                    activeType;

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
            else
            {
                FatalErrorIn
                (
                    "void extendedBlockLduMatrix<tensor>::mapOffDiagCoeffs\n"
                    "(\n"
                    "    const BlockLduMatrix<tensor>& blockLdum\n"
                    ")"
                )   << "Problem between ordinary block matrix and extended"
                    << " block matrix upper coeffs type morphing."
                    << abort(FatalError);
            }
        }
        else
        {
            // Allocate extended lower only
            extendedLowerPtr_ = new TypeCoeffField
            (
                extLduAddr_.extendedLowerAddr().size()
            );
            TypeCoeffField& extLower = *extendedLowerPtr_;

            // Get lower coeffs from underlying lduMatrix
            const TypeCoeffField& lower = blockLdum.lower();

            if (lower.activeType() == blockCoeffBase::SCALAR)
            {
                // Helper type definition
                typedef CoeffField<tensor>::scalarTypeField
                    activeType;

                // Get references to fields
                const activeType& activeLower = lower.asScalar();
                activeType& activeExtLower = extLower.asScalar();

                // Copy non-zero coeffs from basic lduMatrix into corresponding
                // positions
                forAll (lower, faceI)
                {
                    activeExtLower[faceMap[faceI]] = activeLower[faceI];
                }
            }
            else if (lower.activeType() == blockCoeffBase::LINEAR)
            {
                // Helper type definition
                typedef CoeffField<tensor>::linearTypeField
                    activeType;

                // Get references to fields
                const activeType& activeLower = lower.asLinear();
                activeType& activeExtLower = extLower.asLinear();

                // Copy non-zero coeffs from basic lduMatrix into corresponding
                // positions
                forAll (lower, faceI)
                {
                    activeExtLower[faceMap[faceI]] = activeLower[faceI];
                }
            }
            else
            {
                FatalErrorIn
                (
                    "void extendedBlockLduMatrix<tensor>::mapOffDiagCoeffs\n"
                    "(\n"
                    "    const BlockLduMatrix<tensor>& blockLdum\n"
                    ")"
                )   << "Problem between ordinary block matrix and extended"
                    << " block matrix lower coeffs type morphing."
                    << abort(FatalError);
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
            typedef CoeffField<tensor>::scalarTypeField activeType;

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
            typedef CoeffField<tensor>::linearTypeField activeType;

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
        else
        {
            FatalErrorIn
            (
                "void extendedBlockLduMatrix<tensor>::mapOffDiagCoeffs\n"
                "(\n"
                "    const BlockLduMatrix<tensor>& blockLdum\n"
                ")"
            )   << "Problem between ordinary block matrix and extended"
                << " block matrix upper/lower coeffs type morphing."
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
