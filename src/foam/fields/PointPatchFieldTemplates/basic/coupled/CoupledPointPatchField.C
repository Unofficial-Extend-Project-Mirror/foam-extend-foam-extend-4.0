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

#include "CoupledPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class CoupledPointPatch,
    template<class> class MatrixType,
    class Type
>
CoupledPointPatchField
<PatchField, Mesh, PointPatch, CoupledPointPatch, MatrixType, Type>::
CoupledPointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    PatchField<Type>(p, iF)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class CoupledPointPatch,
    template<class> class MatrixType,
    class Type
>
CoupledPointPatchField
<PatchField, Mesh, PointPatch, CoupledPointPatch, MatrixType, Type>::
CoupledPointPatchField
(
    const CoupledPointPatchField
        <PatchField, Mesh, PointPatch, CoupledPointPatch, MatrixType, Type>&,
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const PointPatchFieldMapper&
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    PatchField<Type>(p, iF)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class CoupledPointPatch,
    template<class> class MatrixType,
    class Type
>
CoupledPointPatchField
<PatchField, Mesh, PointPatch, CoupledPointPatch, MatrixType, Type>::
CoupledPointPatchField
(
    const CoupledPointPatchField
    <PatchField, Mesh, PointPatch, CoupledPointPatch, MatrixType, Type>& ptf
)
:
    lduInterfaceField(refCast<const lduInterface>(ptf.patch())),
    PatchField<Type>(ptf)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class CoupledPointPatch,
    template<class> class MatrixType,
    class Type
>
CoupledPointPatchField
<PatchField, Mesh, PointPatch, CoupledPointPatch, MatrixType, Type>::
CoupledPointPatchField
(
    const CoupledPointPatchField
    <PatchField, Mesh, PointPatch, CoupledPointPatch, MatrixType, Type>& ptf,
    const DimensionedField<Type, Mesh>& iF
)
:
    lduInterfaceField(refCast<const lduInterface>(ptf.patch())),
    PatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
