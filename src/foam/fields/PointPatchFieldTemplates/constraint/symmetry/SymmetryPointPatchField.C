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

#include "SymmetryPointPatchField.H"
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
    class SymmetryPointPatch,
    template<class> class MatrixType,
    class Type
>
SymmetryPointPatchField
<PatchField, Mesh, PointPatch, SymmetryPointPatch, MatrixType, Type>::
SymmetryPointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF
)
:
    BasicSymmetryPointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(p, iF)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class SymmetryPointPatch,
    template<class> class MatrixType,
    class Type
>
SymmetryPointPatchField
<PatchField, Mesh, PointPatch, SymmetryPointPatch, MatrixType, Type>::
SymmetryPointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const dictionary& dict
)
:
    BasicSymmetryPointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(p, iF)
{
    if (!isType<SymmetryPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "SymmetryPointPatchField"
            "<PatchField, Mesh, PointPatch, SymmetryPointPatch, "
            "MatrixType, Type>::SymmetryPointPatchField\n"
            "(\n"
            "    const PointPatch& p,\n"
            "    const DimensionedField<Type, Mesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not symmetry type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class SymmetryPointPatch,
    template<class> class MatrixType,
    class Type
>
SymmetryPointPatchField
<PatchField, Mesh, PointPatch, SymmetryPointPatch, MatrixType, Type>::
SymmetryPointPatchField
(
    const SymmetryPointPatchField
        <PatchField, Mesh, PointPatch, SymmetryPointPatch, MatrixType, Type>&,
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const PointPatchFieldMapper&
)
:
    BasicSymmetryPointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(p, iF)
{
    if (!isType<SymmetryPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "SymmetryPointPatchField"
            "<PatchField, Mesh, PointPatch, SymmetryPointPatch, "
            "MatrixType, Type>::SymmetryPointPatchField\n"
            "(\n"
            "    const SymmetryPointPatchField"
                 "<PatchField, Mesh, PointPatch, SymmetryPointPatch, "
            "MatrixType, Type>& ptf,\n"
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
    class SymmetryPointPatch,
    template<class> class MatrixType,
    class Type
>
SymmetryPointPatchField
<PatchField, Mesh, PointPatch, SymmetryPointPatch, MatrixType, Type>::
SymmetryPointPatchField
(
    const SymmetryPointPatchField
    <PatchField, Mesh, PointPatch, SymmetryPointPatch, MatrixType, Type>& ptf
)
:
    BasicSymmetryPointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(ptf)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class SymmetryPointPatch,
    template<class> class MatrixType,
    class Type
>
SymmetryPointPatchField
<PatchField, Mesh, PointPatch, SymmetryPointPatch, MatrixType, Type>::
SymmetryPointPatchField
(
    const SymmetryPointPatchField
    <PatchField, Mesh, PointPatch, SymmetryPointPatch, MatrixType, Type>& ptf,
    const DimensionedField<Type, Mesh>& iF
)
:
    BasicSymmetryPointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
