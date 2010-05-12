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

#include "error.H"

#include "TimeVaryingUniformFixedValuePointPatchField.H"
#include "Time.H"
#include "IFstream.H"

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
TimeVaryingUniformFixedValuePointPatchField
<PatchField, Mesh, PointPatch, MatrixType, Type>::
TimeVaryingUniformFixedValuePointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF
)
:
    FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(p, iF)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
TimeVaryingUniformFixedValuePointPatchField
<PatchField, Mesh, PointPatch, MatrixType, Type>::
TimeVaryingUniformFixedValuePointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const dictionary& dict
)
:
    FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(p, iF),
    timeSeries_(dict)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
TimeVaryingUniformFixedValuePointPatchField
<PatchField, Mesh, PointPatch, MatrixType, Type>::
TimeVaryingUniformFixedValuePointPatchField
(
    const
    TimeVaryingUniformFixedValuePointPatchField
    <PatchField, Mesh, PointPatch, MatrixType, Type>& ptf,
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(p, iF),
    timeSeries_(ptf.timeSeries_)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
TimeVaryingUniformFixedValuePointPatchField
<PatchField, Mesh, PointPatch, MatrixType, Type>::
TimeVaryingUniformFixedValuePointPatchField
(
    const TimeVaryingUniformFixedValuePointPatchField
    <PatchField, Mesh, PointPatch, MatrixType, Type>& ptf
)
:
    FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(ptf),
    timeSeries_(ptf.timeSeries_)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
TimeVaryingUniformFixedValuePointPatchField
<PatchField, Mesh, PointPatch, MatrixType, Type>::
TimeVaryingUniformFixedValuePointPatchField
(
    const TimeVaryingUniformFixedValuePointPatchField
    <PatchField, Mesh, PointPatch, MatrixType, Type>& ptf,
    const DimensionedField<Type, Mesh>& iF
)
:
    FixedValuePointPatchField
        <PatchField, Mesh, PointPatch, MatrixType, Type>(ptf, iF),
    timeSeries_(ptf.timeSeries_)
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
TimeVaryingUniformFixedValuePointPatchField
<PatchField, Mesh, PointPatch, MatrixType, Type>::
initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    Field<Type>& patchField = *this;

    patchField = timeSeries_(this->db().time().timeOutputValue());

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
void TimeVaryingUniformFixedValuePointPatchField
<PatchField, Mesh, PointPatch, MatrixType, Type>::write(Ostream& os) const
{
    PatchField<Type>::write(os);
    timeSeries_.write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
