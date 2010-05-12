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
