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

#include "immersedBoundaryKqRWallFunctionFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "standAlonePatch.H"
#include "surfaceWriter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void immersedBoundaryKqRWallFunctionFvPatchField<Type>::checkType()
{
    if (!this->patch().isWall())
    {
        FatalErrorIn("immersedBoundaryKqRWallFunctionFvPatchField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << this->patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << this->patch().type()
            << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
immersedBoundaryKqRWallFunctionFvPatchField<Type>::
immersedBoundaryKqRWallFunctionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    kqRWallFunctionFvPatchField<Type>(p, iF),
    immersedBoundaryFieldBase<Type>(p, true, pTraits<Type>::one*SMALL)
{
    this->checkType();
}


template<class Type>
immersedBoundaryKqRWallFunctionFvPatchField<Type>::
immersedBoundaryKqRWallFunctionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    kqRWallFunctionFvPatchField<Type>(p, iF),   // Do not read mixed data
    immersedBoundaryFieldBase<Type>
    (
        p,
        Switch(dict.lookup("setDeadValue")),
        pTraits<Type>(dict.lookup("deadValue"))
    )
{
    this->readPatchType(dict);
    this->checkType();

    Field<Type>::operator=(this->patchInternalField());
}


template<class Type>
immersedBoundaryKqRWallFunctionFvPatchField<Type>::
immersedBoundaryKqRWallFunctionFvPatchField
(
    const immersedBoundaryKqRWallFunctionFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    kqRWallFunctionFvPatchField<Type>(p, iF), // Do not map mixed data.  Set patchType later
    immersedBoundaryFieldBase<Type>
    (
        p,
        ptf.setDeadValue(),
        ptf.deadValue()
    )
{
    this->setPatchType(ptf);
    this->checkType();

    // cannot be used.  Initialise the value to avoid errors
    // HJ, 1/Dec/2017
    Field<Type>::operator=(pTraits<Type>::zero);
}


template<class Type>
immersedBoundaryKqRWallFunctionFvPatchField<Type>::
immersedBoundaryKqRWallFunctionFvPatchField
(
    const immersedBoundaryKqRWallFunctionFvPatchField& ptf
)
:
    kqRWallFunctionFvPatchField<Type>(ptf),
    immersedBoundaryFieldBase<Type>
    (
        ptf.ibPatch(),
        ptf.setDeadValue(),
        ptf.deadValue()
    )
{
    this->setPatchType(ptf);
    this->checkType();
}


template<class Type>
immersedBoundaryKqRWallFunctionFvPatchField<Type>::
immersedBoundaryKqRWallFunctionFvPatchField
(
    const immersedBoundaryKqRWallFunctionFvPatchField& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    kqRWallFunctionFvPatchField<Type>(ptf, iF),
    immersedBoundaryFieldBase<Type>
    (
        ptf.ibPatch(),
        ptf.setDeadValue(),
        ptf.deadValue()
    )
{
    this->setPatchType(ptf);
    this->checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void immersedBoundaryKqRWallFunctionFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper&
)
{
    Field<Type>::operator=(this->patchInternalField());
}


template<class Type>
void immersedBoundaryKqRWallFunctionFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList&
)
{}


template<class Type>
void immersedBoundaryKqRWallFunctionFvPatchField<Type>::updateOnMotion()
{
    if (this->size() != this->ibPatch().size())
    {
        // Use internal values, resizing the file if needed
        Field<Type>::operator=(this->patchInternalField());
    }
}


template<class Type>
void immersedBoundaryKqRWallFunctionFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Resize on evaluation if needed
    if (this->size() != this->ibPatch().size())
    {
        // Use internal values, resizing the file if needed
        Field<Type>::operator=(this->patchInternalField());
    }

    // Get non-constant reference to internal field
    Field<Type>& intField = const_cast<Field<Type>&>(this->internalField());

    // Set dead value
    this->setDeadValues(intField);

    kqRWallFunctionFvPatchField<Type>::evaluate(commsType);
}


template<class Type>
void immersedBoundaryKqRWallFunctionFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix
)
{
    this->setDeadValues(matrix);
}


template<class Type>
void immersedBoundaryKqRWallFunctionFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeDeadData(os);

    // The value entry needs to be written with zero size
    Field<Type>::null().writeEntry("value", os);
    // this->writeEntry("value", os);

    this->writeField(*this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
