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

\*---------------------------------------------------------------------------*/

#include "fixedValueIbFvPatchField.H"
#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void fixedValueIbFvPatchField<Type>::updateIbValues()
{
    // Interpolate the values from tri surface using nearest triangle
    const labelList& nt = this->ibPatch().ibPolyPatch().nearestTri();

    Field<Type>::operator=(Field<Type>(triValue_, nt));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
fixedValueIbFvPatchField<Type>::fixedValueIbFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    immersedBoundaryFieldBase<Type>(p, false, pTraits<Type>::zero),
    triValue_(this->ibPatch().ibMesh().size(), pTraits<Type>::zero)
{}


template<class Type>
fixedValueIbFvPatchField<Type>::fixedValueIbFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),   // Do not read mixed data
    immersedBoundaryFieldBase<Type>
    (
        p,
        Switch(dict.lookup("setDeadValue")),
        pTraits<Type>(dict.lookup("deadValue"))
    ),
    triValue_("triValue", dict, this->ibPatch().ibMesh().size())
{
    // Since patch does not read a dictionary, the patch type needs to be read
    // manually.  HJ, 6/Sep/2018
    this->readPatchType(dict);

    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    // Re-interpolate the data related to immersed boundary
    this->updateIbValues();

    fixedValueFvPatchField<Type>::evaluate();
}


template<class Type>
fixedValueIbFvPatchField<Type>::fixedValueIbFvPatchField
(
    const fixedValueIbFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(p, iF),  // Do not map fixed data
    immersedBoundaryFieldBase<Type>
    (
        p,
        ptf.setDeadValue(),
        ptf.deadValue()
    ),
    triValue_(ptf.triValue())
{
    // Note: NO MAPPING.  Fields are created on the immersed boundary
    // HJ, 12/Apr/2012
    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalErrorInFunction
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    // Re-interpolate the data related to immersed boundary
    this->updateIbValues();
    this->setPatchType(ptf);

    // On creation of the mapped field, the internal field is dummy and
    // cannot be used.  Initialise the value to avoid errors
    // HJ, 1/Dec/2017
    Field<Type>::operator=(pTraits<Type>::zero);
}


template<class Type>
fixedValueIbFvPatchField<Type>::fixedValueIbFvPatchField
(
    const fixedValueIbFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    immersedBoundaryFieldBase<Type>
    (
        ptf.ibPatch(),
        ptf.setDeadValue(),
        ptf.deadValue()
    ),
    triValue_(ptf.triValue())
{
    this->setPatchType(ptf);
}


template<class Type>
fixedValueIbFvPatchField<Type>::fixedValueIbFvPatchField
(
    const fixedValueIbFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    immersedBoundaryFieldBase<Type>
    (
        ptf.ibPatch(),
        ptf.setDeadValue(),
        ptf.deadValue()
    ),
    triValue_(ptf.triValue())
{
    this->setPatchType(ptf);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void fixedValueIbFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    // Base fields do not map: re-interpolate them from tri data
    this->updateIbValues();
}


template<class Type>
void fixedValueIbFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList&
)
{
    // Base fields do not rmap: re-interpolate them from tri data

    const fixedValueIbFvPatchField<Type>& mptf =
        refCast<const fixedValueIbFvPatchField<Type> >(ptf);

    // Set rmap tri data
    triValue_ = mptf.triValue_;

    this->updateIbValues();
}


template<class Type>
void fixedValueIbFvPatchField<Type>::updateOnMotion()
{
    if (this->size() != this->ibPatch().size())
    {
        this->updateIbValues();
    }
}


template<class Type>
void fixedValueIbFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    this->updateIbValues();

    // Get non-constant reference to internal field
    Field<Type>& intField = const_cast<Field<Type>&>(this->internalField());

    // Set dead value
    this->setDeadValues(intField);

    // Evaluate fixed value condition
    fixedValueFvPatchField<Type>::evaluate();
}


template<class Type>
void fixedValueIbFvPatchField<Type>::write(Ostream& os) const
{
    // Resolve post-processing issues.  HJ, 1/Dec/2017
    fvPatchField<Type>::write(os);
    triValue_.writeEntry("triValue", os);
    immersedBoundaryFieldBase<Type>::writeDeadData(os);

    // The value entry needs to be written with zero size
    Field<Type>::null().writeEntry("value", os);
    // this->writeEntry("value", os);

    this->writeField(*this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
