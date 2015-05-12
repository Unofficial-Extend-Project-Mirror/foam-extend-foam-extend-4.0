/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

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
