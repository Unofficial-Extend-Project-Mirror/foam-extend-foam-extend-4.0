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

Description

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "UniformFixedValuePointPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
UniformFixedValuePointPatchField
<PatchField, Mesh, PointPatch, MatrixType, Type>::
UniformFixedValuePointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF
)
:
    FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(p, iF),
    uniformValue_(pTraits<Type>::zero)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
UniformFixedValuePointPatchField
<PatchField, Mesh, PointPatch, MatrixType, Type>::
UniformFixedValuePointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const dictionary& dict
)
:
    FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(p, iF),
    uniformValue_(pTraits<Type>(dict.lookup("uniformValue")))
{
    FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>::operator==
        (uniformValue_);
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
UniformFixedValuePointPatchField
<PatchField, Mesh, PointPatch, MatrixType, Type>::
UniformFixedValuePointPatchField
(
    const
    UniformFixedValuePointPatchField
    <PatchField, Mesh, PointPatch, MatrixType, Type>& ptf,
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(p, iF),
    uniformValue_(ptf.uniformValue_)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
UniformFixedValuePointPatchField
<PatchField, Mesh, PointPatch, MatrixType, Type>::
UniformFixedValuePointPatchField
(
    const UniformFixedValuePointPatchField
    <PatchField, Mesh, PointPatch, MatrixType, Type>& ptf
)
:
    FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(ptf),
    uniformValue_(ptf.uniformValue_)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
UniformFixedValuePointPatchField
<PatchField, Mesh, PointPatch, MatrixType, Type>::
UniformFixedValuePointPatchField
(
    const UniformFixedValuePointPatchField
    <PatchField, Mesh, PointPatch, MatrixType, Type>& ptf,
    const DimensionedField<Type, Mesh>& iF
)
:
    FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(ptf, iF),
    uniformValue_(ptf.uniformValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
void
UniformFixedValuePointPatchField
<PatchField, Mesh, PointPatch, MatrixType, Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    Field<Type>& patchField = *this;
    patchField = uniformValue_;

    FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>::
        initEvaluate(commsType);
}


// Write
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
void
UniformFixedValuePointPatchField
<PatchField, Mesh, PointPatch, MatrixType, Type>::write(Ostream& os) const
{
    PatchField<Type>::write(os);
    os.writeKeyword("uniformValue")
        << uniformValue_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
