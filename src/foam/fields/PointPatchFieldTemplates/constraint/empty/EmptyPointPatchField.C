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

#include "EmptyPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class EmptyPointPatch,
    template<class> class MatrixType,
    class Type
>
EmptyPointPatchField
<PatchField, Mesh, PointPatch, EmptyPointPatch, MatrixType, Type>::
EmptyPointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF
)
:
    PatchField<Type>(p, iF)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class EmptyPointPatch,
    template<class> class MatrixType,
    class Type
>
EmptyPointPatchField
<PatchField, Mesh, PointPatch, EmptyPointPatch, MatrixType, Type>::
EmptyPointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const dictionary& dict
)
:
    PatchField<Type>(p, iF)
{
    if (!isType<EmptyPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "EmptyPointPatchField"
            "<PatchField, Mesh, PointPatch, EmptyPointPatch, "
            "MatrixType, Type>::EmptyPointPatchField\n"
            "(\n"
            "    const PointPatch& p,\n"
            "    const DimensionedField<Type, Mesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not empty type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class EmptyPointPatch,
    template<class> class MatrixType,
    class Type
>
EmptyPointPatchField
<PatchField, Mesh, PointPatch, EmptyPointPatch, MatrixType, Type>::
EmptyPointPatchField
(
    const EmptyPointPatchField
    <PatchField, Mesh, PointPatch, EmptyPointPatch, MatrixType, Type>&,
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const PointPatchFieldMapper&
)
:
    PatchField<Type>(p, iF)
{
    if (!isType<EmptyPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "EmptyPointPatchField"
            "<PatchField, Mesh, PointPatch, EmptyPointPatch, "
            "MatrixType, Type>::EmptyPointPatchField\n"
            "(\n"
            "    const EmptyPointPatchField<PatchField, Mesh, PointPatch, "
                 "EmptyPointPatch, MatrixType, Type>& ptf,\n"
            "    const PointPatch& p,\n"
            "    const DimensionedField<Type, Mesh>& iF,\n"
            "    const PointPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class EmptyPointPatch,
    template<class> class MatrixType,
    class Type
>
EmptyPointPatchField
<PatchField, Mesh, PointPatch, EmptyPointPatch, MatrixType, Type>::
EmptyPointPatchField
(
    const EmptyPointPatchField
    <PatchField, Mesh, PointPatch, EmptyPointPatch, MatrixType, Type>& ptf
)
:
    PatchField<Type>(ptf)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class EmptyPointPatch,
    template<class> class MatrixType,
    class Type
>
EmptyPointPatchField
<PatchField, Mesh, PointPatch, EmptyPointPatch, MatrixType, Type>::
EmptyPointPatchField
(
    const EmptyPointPatchField
    <PatchField, Mesh, PointPatch, EmptyPointPatch, MatrixType, Type>& ptf,
    const DimensionedField<Type, Mesh>& iF
)
:
    PatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
