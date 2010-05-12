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

Description

\*---------------------------------------------------------------------------*/

#include "CalculatedPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
const word&
PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>
::calculatedType()
{
    return CalculatedPointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>::typeName;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
CalculatedPointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
CalculatedPointPatchField
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
    template<class> class MatrixType,
    class Type
>
CalculatedPointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
CalculatedPointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const dictionary&
)
:
    PatchField<Type>(p, iF)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
CalculatedPointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
CalculatedPointPatchField
(
    const CalculatedPointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>&,
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const PointPatchFieldMapper&
)
:
    PatchField<Type>(p, iF)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
CalculatedPointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
CalculatedPointPatchField
(
    const CalculatedPointPatchField
    <PatchField, Mesh, PointPatch, MatrixType, Type>& ptf
)
:
    PatchField<Type>(ptf)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
CalculatedPointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
CalculatedPointPatchField
(
    const CalculatedPointPatchField
    <PatchField, Mesh, PointPatch, MatrixType, Type>& ptf,
    const DimensionedField<Type, Mesh>& iF
)
:
    PatchField<Type>(ptf, iF)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
template<class Type2>
autoPtr<PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type> >
PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>
::NewCalculatedType
(
    const PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type2>& pf
)
{
    typename PointPatchConstructorTable::iterator patchTypeCstrIter =
        PointPatchConstructorTablePtr_->find(pf.patch().type());

    if (patchTypeCstrIter != PointPatchConstructorTablePtr_->end())
    {
        return autoPtr<PatchField<Type> >
        (
            patchTypeCstrIter()
            (
                pf.patch(),
                DimensionedField<Type, Mesh>::null()
            )
        );
    }
    else
    {
        return
            autoPtr
            <
                PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>
            >
            (
                new CalculatedPointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>
                (
                    pf.patch(),
                    DimensionedField<Type, Mesh>::null()
                )
            );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
