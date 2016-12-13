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

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::Keyed<T>::Keyed()
:
    key_(-1)
{}


template<class T>
inline Foam::Keyed<T>::Keyed(const T& item, const label key)
:
    T(item),
    key_(key)
{}


template<class T>
inline Foam::Keyed<T>::Keyed(const Xfer<T>& item, const label key)
:
    T(item),
    key_(key)
{}


template<class T>
inline Foam::Keyed<T>::Keyed(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::Keyed<T>::key() const
{
    return key_;
}

template<class T>
inline Foam::label& Foam::Keyed<T>::key()
{
    return key_;
}


template<class T>
inline Foam::List<Foam::Keyed<T> >
Foam::Keyed<T>::createList(const List<T>& lst, const label key)
{
    List<Keyed<T> > newList(lst.size());

    forAll (lst, elemI)
    {
        newList[elemI] = Keyed(lst[elemI], key);
    }
    return newList;
}


template<class T>
inline Foam::List<Foam::Keyed<T> >
Foam::Keyed<T>::createList(const List<T>& lst, const List<label>& keys)
{
    if (lst.size() != keys.size())
    {
        FatalErrorIn
        (
            "Foam::Keyed<T>::createList(const List<T>&, const List<label>&)"
        )
            << "size mismatch adding keys to a list:" << nl
            << "List has size " << lst.size()
            << " and keys has size " << keys.size() << nl
            << abort(FatalError);
    }

    List<Keyed<T> > newList(lst.size());

    forAll (lst, elemI)
    {
        newList[elemI] = Keyed(lst[elemI], keys[elemI]);
    }
    return newList;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class T>
inline Foam::Istream& Foam::operator>>(Istream& is, Keyed<T>& item)
{
    // Read beginning of Keyed item/key pair
    is.readBegin("Keyed<T>");

    is >> static_cast<T&>(item);
    is >> item.key_;

    // Read end of Keyed item/key pair
    is.readEnd("Keyed<T>");

    // Check state of Ostream
    is.check("Istream& operator>>(Istream&, Keyed&)");

    return is;
}


template<class T>
inline Foam::Ostream& Foam::operator<<(Ostream& os, const Keyed<T>& item)
{
    os  << token::BEGIN_LIST
        << static_cast<const T&>(item)
        << token::SPACE << item.key_
        << token::END_LIST;

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
