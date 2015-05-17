/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "immersedBoundaryWallFunctionFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void immersedBoundaryWallFunctionFvPatchField<Type>::motionUpdate() const
{
    wallValue_.clear();
    wallMask_.clear();

    immersedBoundaryFvPatchField<Type>::motionUpdate();
}


template<class Type>
void immersedBoundaryWallFunctionFvPatchField<Type>::setIbCellValues
(
    const Field<Type>& ibcValues
) const
{
    const labelList& ibc = this->ibPatch().ibCells();

    if (ibcValues.size() != ibc.size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "void immersedBoundaryWallFunctionFvPatchField<Type>::"
            "setIbCellValues\n"
            "(\n"
            "    const Field<Type>& ibcValues\n"
            ") const"
        )   << "Size of ibcValues field not equal to the number of IB cells."
            << nl << "ibcValues: " << ibcValues.size()
            << " ibc: " << ibc.size()
            << abort(FatalError);
    }

    // Get non-const access to internal field
    Field<Type>& psiI = const_cast<Field<Type>&>(this->internalField());

    if (wallValue_.empty() || wallMask_.empty())
    {
        immersedBoundaryFvPatchField<Type>::setIbCellValues(ibcValues);
    }
    else
    {
        forAll (ibcValues, cellI)
        {
            // If mask is set use the wall value, otherwise use the
            // fitted value
            if (wallMask_[cellI])
            {
                psiI[ibc[cellI]] = wallValue_[cellI];
            }
            else
            {
                psiI[ibc[cellI]] = ibcValues[cellI];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedBoundaryFvPatchField<Type>(p, iF),
    wallValue_(),
    wallMask_()
{}


template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    immersedBoundaryFvPatchField<Type>(p, iF, dict),
    wallValue_(),
    wallMask_()
{}


template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const immersedBoundaryWallFunctionFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedBoundaryFvPatchField<Type>(ptf, p, iF, mapper),
    wallValue_(),
    wallMask_()
{}


template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const immersedBoundaryWallFunctionFvPatchField& tkqrwfpf
)
:
    immersedBoundaryFvPatchField<Type>(tkqrwfpf),
    wallValue_(),
    wallMask_()
{}


template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const immersedBoundaryWallFunctionFvPatchField& tkqrwfpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedBoundaryFvPatchField<Type>(tkqrwfpf, iF),
    wallValue_(),
    wallMask_()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
