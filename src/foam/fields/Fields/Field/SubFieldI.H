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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::SubField<Type>::SubField
(
    const SubList<Type>& list
)
:
    SubList<Type>(list)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const UList<Type>& list
)
:
    SubList<Type>(list, list.size())
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const UList<Type>& list,
    const label subSize
)
:
    SubList<Type>(list, subSize)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const UList<Type>& list,
    const label subSize,
    const label startIndex
)
:
    SubList<Type>(list, subSize, startIndex)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const SubField<Type>& sfield
)
:
    refCount(),
    SubList<Type>(sfield)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline const Foam::SubField<Type>& Foam::SubField<Type>::null()
{
    // Dangerous, but no null constructor.  Reconsider.  HJ, 10/Apr/2013
    SubField<Type>* nullPtr = reinterpret_cast<SubField<Type>*>(NULL);
    return *nullPtr;
}


template<class Type>
inline Foam::tmp<Foam::Field<typename Foam::SubField<Type>::cmptType> >
Foam::SubField<Type>::component
(
    const direction d
) const
{
    return (reinterpret_cast<const Field<Type>&>(*this)).component(d);
}


template<class Type>
inline Foam::tmp<Foam::Field<Type> > Foam::SubField<Type>::T() const
{
    return (reinterpret_cast<const Field<Type>&>(*this)).T();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline void Foam::SubField<Type>::operator=(const SubField<Type>& rhs)
{
    UList<Type>::operator=(rhs);
}


template<class Type>
inline Foam::SubField<Type>::operator const Foam::Field<Type>&() const
{
    return *reinterpret_cast< const Field<Type>* >(this);
}


// ************************************************************************* //
