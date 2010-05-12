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
autoPtr<PatchField<Type> >
PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::New
(
    const word& patchFieldType,
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF
)
{
    if (debug)
    {
        InfoIn
        (
            "PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::"
            "New(const word&, const PointPatch&, "
            "const DimensionedField<Type, Mesh>&)"
        )   << "constructing PointPatchField<PatchField, PointPatch, "
            << "MatrixType, Type>"
            << endl;
    }

    typename PointPatchConstructorTable::iterator cstrIter =
        PointPatchConstructorTablePtr_->find(patchFieldType);

    if (cstrIter == PointPatchConstructorTablePtr_->end())
    {
        cstrIter = PointPatchConstructorTablePtr_->find("default");

        if (cstrIter == PointPatchConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "PointPatchField<PatchField, Mesh, PointPatch, "
                "MatrixType, Type>::"
                "New(const word&, const PointPatch&, "
                "const DimensionedField<Type, Mesh>&)"
            )   << "Unknown patchTypefield type "
                << patchFieldType
                << endl << endl
                << "Valid patchField types are :" << endl
                << PointPatchConstructorTablePtr_->toc()
                << exit(FatalError);
        }
    }

    typename PointPatchConstructorTable::iterator patchTypeCstrIter =
        PointPatchConstructorTablePtr_->find(p.type());

    if (patchTypeCstrIter != PointPatchConstructorTablePtr_->end())
    {
        return autoPtr<PatchField<Type> >(patchTypeCstrIter()(p, iF));
    }
    else
    {
        return autoPtr<PatchField<Type> >(cstrIter()(p, iF));
    }
}


// Return a pointer to a new patch created on freestore from
// a given PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>
// mapped onto a new patch
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
autoPtr<PatchField<Type> >
PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::New
(
    const PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>& ptf,
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const PointPatchFieldMapper& pfMapper
)
{
    if (debug)
    {
        InfoIn
        (
            "PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::"
            "New(const PointPatchField<PatchField, PointPatch, MatrixType, "
            "Type>&, const PointPatch&, const DimensionedField<Type, Mesh>&, "
            "const PointPatchFieldMapper&)"
        )   << "constructing PointPatchField<PatchField, PointPatch, "
            << "MatrixType, Type>"
            << endl;
    }

    typename patchMapperConstructorTable::iterator cstrIter =
        patchMapperConstructorTablePtr_->find(ptf.type());

    if (cstrIter == patchMapperConstructorTablePtr_->end())
    {
        cstrIter = patchMapperConstructorTablePtr_->find("default");

        if (cstrIter == patchMapperConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "PointPatchField<PatchField, Mesh, PointPatch, "
                "MatrixType, Type>::"
                "New(const PointPatchField<PatchField, PointPatch, "
                "MatrixType, Type>&, const PointPatch&, "
                "const DimensionedField<Type, Mesh>&, "
                "const PointPatchFieldMapper&)"
            )   << "unknown patchTypefield type "
                << ptf.type() << endl << endl
                << "Valid patchField types are :" << endl
                << patchMapperConstructorTablePtr_->toc()
                << exit(FatalError);
        }
    }

    typename patchMapperConstructorTable::iterator
        patchTypeCstrIter = patchMapperConstructorTablePtr_->find(p.type());

    if (patchTypeCstrIter != patchMapperConstructorTablePtr_->end())
    {
        return autoPtr<PatchField<Type> >
        (
            patchTypeCstrIter()(ptf, p, iF, pfMapper)
        );
    }
    else
    {
        return autoPtr<PatchField<Type> >(cstrIter()(ptf, p, iF, pfMapper));
    }
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
autoPtr<PatchField<Type> >
PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::New
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const dictionary& dict
)
{
    if (debug)
    {
        InfoIn
        (
            "PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::"
            "New(const PointPatch&, const DimensionedField<Type, Mesh>&, "
            "const dictionary&)"
        )   << "constructing PointPatchField<PatchField, PointPatch, "
            << "MatrixType, Type>"
            << endl;
    }

    word patchFieldType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter
        = dictionaryConstructorTablePtr_->find(patchFieldType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        cstrIter = dictionaryConstructorTablePtr_->find("default");

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "PointPatchField<PatchField, Mesh, PointPatch, "
                "MatrixType, Type>::"
                "New(const PointPatch&, const DimensionedField<Type, Mesh>&, "
                "const dictionary&)",
                dict
            )   << "Unknown patchField type " << patchFieldType
                << " for patch type " << p.type() << endl << endl
                << "Valid patchField types are :" << endl
                << dictionaryConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }
    }

    typename dictionaryConstructorTable::iterator patchTypeCstrIter
        = dictionaryConstructorTablePtr_->find(p.type());

    if
    (
        patchTypeCstrIter != dictionaryConstructorTablePtr_->end()
     && patchTypeCstrIter() != cstrIter()
    )
    {
        FatalIOErrorIn
        (
            "PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>"
            "const PointPatch&, const DimensionedField<Type, Mesh>&, "
            "const dictionary&)",
            dict
        )   << "inconsistent patch and patchField types for \n"
            << "    patch type " << p.type()
            << " and patchField type " << patchFieldType
            << exit(FatalIOError);
    }

    return autoPtr<PatchField<Type> >(cstrIter()(p, iF, dict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
