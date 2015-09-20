/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faePatchField<Type> > faePatchField<Type>::New
(
    const word& patchFieldType,
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF
)
{
    if (debug)
    {
        Info<< "faePatchField<Type>::New(const word&, const faPatch&, "
               "const DimensionedField<Type, edgeMesh>&) : "
               "constructing faePatchField<Type>"
            << endl;
    }

    typename patchConstructorTable::iterator cstrIter =
        patchConstructorTablePtr_->find(patchFieldType);

    if (cstrIter == patchConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "faePatchField<Type>::New(const word&, const faPatch&, "
            "const DimensionedField<Type, edgeMesh>&)"
        )   << "Unknown patchTypefield type " << patchFieldType
            << endl << endl
            << "Valid patchField types are :" << endl
            << patchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    typename patchConstructorTable::iterator patchTypeCstrIter =
        patchConstructorTablePtr_->find(p.type());

    if (patchTypeCstrIter != patchConstructorTablePtr_->end())
    {
        return patchTypeCstrIter()(p, iF);
    }
    else
    {
        return cstrIter()(p, iF);
    }
}


template<class Type>
tmp<faePatchField<Type> > faePatchField<Type>::New
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const dictionary& dict
)
{
    if (debug)
    {
        Info<< "faePatchField<Type>::New(const faPatch&, "
               "const DimensionedField<Type, edgeMesh>&, "
               "const dictionary&) : "
               "constructing faePatchField<Type>"
            << endl;
    }

    word patchFieldType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter
        = dictionaryConstructorTablePtr_->find(patchFieldType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        if (!disallowDefaultFaePatchField)
        {
            cstrIter = dictionaryConstructorTablePtr_->find("default");
        }

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "faePatchField<Type>::New(const faPatch&, "
                "const DimensionedField<Type, edgeMesh>&, "
                "const dictionary&)",
                dict
            )   << "Unknown patchField type " << patchFieldType
                << " for patch type " << p.type() << endl << endl
                << "Valid patchField types are :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }
    }

    typename dictionaryConstructorTable::iterator patchTypeCstrIter
        = dictionaryConstructorTablePtr_->find(p.type());

    if
    (
        patchTypeCstrIter != dictionaryConstructorTablePtr_->end()
     && *patchTypeCstrIter != *cstrIter
    )
    {
        FatalIOErrorIn
        (
            "faePatchField<Type>const faPatch&, "
            "const DimensionedField<Type, edgeMesh>&, "
            "const dictionary&)",
            dict
        ) << "inconsistent patch and patchField types for \n"
             "    patch type " << p.type()
            << " and patchField type " << patchFieldType
            << exit(FatalIOError);
    }

    return cstrIter()(p, iF, dict);
}


// Return a pointer to a new patch created on freestore from
// a given faePatchField<Type> mapped onto a new patch
template<class Type>
tmp<faePatchField<Type> > faePatchField<Type>::New
(
    const faePatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const faPatchFieldMapper& pfMapper
)
{
    if (debug)
    {
        Info<< "faePatchField<Type>::New(const faePatchField<Type>&,"
               " const faPatch&, const DimensionedField<Type, edgeMesh>&, "
               "const faPatchFieldMapper&) : "
               "constructing faePatchField<Type>"
            << endl;
    }

    typename patchMapperConstructorTable::iterator cstrIter =
        patchMapperConstructorTablePtr_->find(ptf.type());

    if (cstrIter == patchMapperConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "faePatchField<Type>::New(const faePatchField<Type>&, "
            "const faPatch&, const DimensionedField<Type, edgeMesh>&, "
            "const faPatchFieldMapper&)"
        )   << "unknown patchTypefield type " << ptf.type() << endl << endl
            << "Valid patchField types are :" << endl
            << patchMapperConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    typename patchMapperConstructorTable::iterator
        patchTypeCstrIter = patchMapperConstructorTablePtr_->find(p.type());

    if (patchTypeCstrIter != patchMapperConstructorTablePtr_->end())
    {
        return patchTypeCstrIter()(ptf, p, iF, pfMapper);
    }
    else
    {
        return cstrIter()(ptf, p, iF, pfMapper);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
