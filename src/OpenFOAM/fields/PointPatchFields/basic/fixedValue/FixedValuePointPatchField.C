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

#include "FixedValuePointPatchField.H"
#include "boolList.H"
#include "Map.H"

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
FixedValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
FixedValuePointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF
)
:
    ValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>(p, iF)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
FixedValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
FixedValuePointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const dictionary& dict
)
:
    ValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>
        (p, iF, dict)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
FixedValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
FixedValuePointPatchField
(
    const FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>& ptf,
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    ValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>
        (ptf, p, iF, mapper)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
FixedValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
FixedValuePointPatchField
(
    const FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>& ptf
)
:
    ValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>(ptf)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
FixedValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
FixedValuePointPatchField
(
    const FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>& ptf,
    const DimensionedField<Type, Mesh>& iF
)
:
    ValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>
        (ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set boundary condition to matrix
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
void
FixedValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
setBoundaryCondition
(
    Map<typename MatrixType<Type>::ConstraintType>& fix
) const
{
    const Field<Type>& values = *this;

    // get addressing
    const labelList& meshPoints = this->patch().meshPoints();

    forAll (meshPoints, pointI)
    {
        const label curPoint = meshPoints[pointI];

        // create a constraint
        typename MatrixType<Type>::ConstraintType bc
        (
            curPoint,
            values[pointI],
            pTraits<Type>::one
        );

        // If not set add it, otherwise combine with
        // already existing value
        if (!fix.found(curPoint))
        {
            fix.insert(curPoint, bc);
        }
        else
        {
            fix[curPoint].combine(bc);
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Force an assignment
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
void
FixedValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>
::operator==
(
    const ValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>& ptf
)
{
    Field<Type>::operator=(ptf);

    // Insert the result into the internal field
    initEvaluate();
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
void
FixedValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>
::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);

    // insert the result into the internal field
    initEvaluate();
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
void
FixedValuePointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);

    // Insert the result into the internal field
    initEvaluate();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
