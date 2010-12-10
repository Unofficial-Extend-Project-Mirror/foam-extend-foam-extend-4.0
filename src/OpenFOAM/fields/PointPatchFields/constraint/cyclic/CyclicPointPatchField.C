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

#include "CyclicPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class CyclicPointPatch,
    template<class> class MatrixType,
    class Type
>
CyclicPointPatchField
<PatchField, Mesh, PointPatch, CyclicPointPatch, MatrixType, Type>::
CyclicPointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename CyclicPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(p, iF),
    cyclicPatch_(refCast<const CyclicPointPatch>(p))
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class CyclicPointPatch,
    template<class> class MatrixType,
    class Type
>
CyclicPointPatchField
<PatchField, Mesh, PointPatch, CyclicPointPatch, MatrixType, Type>::
CyclicPointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const dictionary& dict
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh, 
        PointPatch,
        typename CyclicPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(p, iF),
    cyclicPatch_(refCast<const CyclicPointPatch>(p))
{
    if (!isType<CyclicPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "CyclicPointPatchField"
            "<PatchField, Mesh, PointPatch, CyclicPointPatch, "
            "MatrixType, Type>::"
            "CyclicPointPatchField\n"
            "(\n"
            "    const PointPatch& p,\n"
            "    const DimensionedField<Type, Mesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not cyclic type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class CyclicPointPatch,
    template<class> class MatrixType,
    class Type
>
CyclicPointPatchField
<PatchField, Mesh, PointPatch, CyclicPointPatch, MatrixType, Type>::
CyclicPointPatchField
(
    const CyclicPointPatchField
        <PatchField, Mesh, PointPatch, CyclicPointPatch, MatrixType, Type>&,
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const PointPatchFieldMapper&
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename CyclicPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(p, iF),
    cyclicPatch_(refCast<const CyclicPointPatch>(p))
{
    if (!isType<CyclicPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "CyclicPointPatchField"
            "<PatchField, Mesh, PointPatch, CyclicPointPatch, "
            "MatrixType, Type>::CyclicPointPatchField\n"
            "(\n"
            "    const CyclicPointPatchField"
            "    <PatchField, Mesh, PointPatch, CyclicPointPatch, "
            "MatrixType, Type>&"
            "ptf,\n"
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
    class CyclicPointPatch,
    template<class> class MatrixType,
    class Type
>
CyclicPointPatchField
<PatchField, Mesh, PointPatch, CyclicPointPatch, MatrixType, Type>::
CyclicPointPatchField
(
    const CyclicPointPatchField
    <PatchField, Mesh, PointPatch, CyclicPointPatch, MatrixType, Type>& ptf
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename CyclicPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(ptf),
    cyclicPatch_(ptf.cyclicPatch_)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class CyclicPointPatch,
    template<class> class MatrixType,
    class Type
>
CyclicPointPatchField
<PatchField, Mesh, PointPatch, CyclicPointPatch, MatrixType, Type>::
CyclicPointPatchField
(
    const CyclicPointPatchField
    <PatchField, Mesh, PointPatch, CyclicPointPatch, MatrixType, Type>& ptf,
    const DimensionedField<Type, Mesh>& iF
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename CyclicPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(ptf, iF),
    cyclicPatch_(ptf.cyclicPatch_)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
