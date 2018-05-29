/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

inline Foam::PackedBoolList::PackedBoolList()
:
    PackedList<1>()
{}


inline Foam::PackedBoolList::PackedBoolList(const label size)
:
    PackedList<1>(size)
{}


inline Foam::PackedBoolList::PackedBoolList
(
    const label size,
    const bool val
)
:
    PackedList<1>(size, (val ? 1u : 0u))
{}


inline Foam::PackedBoolList::PackedBoolList(const PackedBoolList& lst)
:
    PackedList<1>(lst)
{}


inline Foam::PackedBoolList::PackedBoolList(const PackedList<1>& lst)
:
    PackedList<1>(lst)
{}


inline Foam::PackedBoolList::PackedBoolList(const Xfer<PackedBoolList>& lst)
:
    PackedList<1>()
{
    transfer(lst());
}


inline Foam::PackedBoolList::PackedBoolList(const Xfer<PackedList<1> >& lst)
:
    PackedList<1>(lst)
{}


inline Foam::PackedBoolList::PackedBoolList(const Foam::UList<bool>& lst)
:
    PackedList<1>()
{
    operator=(lst);
}


inline Foam::PackedBoolList::PackedBoolList(const UList<label>& indices)
:
    PackedList<1>(indices.size(), 0u)
{
    set(indices);
}


inline Foam::PackedBoolList::PackedBoolList(const UIndirectList<label>& indices)
:
    PackedList<1>(indices.size(), 0u)
{
    set(indices);
}


inline Foam::autoPtr<Foam::PackedBoolList>
Foam::PackedBoolList::clone() const
{
    return autoPtr<PackedBoolList>(new PackedBoolList(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::PackedBoolList::transfer(PackedBoolList& lst)
{
    PackedList<1>::transfer(static_cast<PackedList<1>&>(lst));
}


inline void Foam::PackedBoolList::transfer(PackedList<1>& lst)
{
    PackedList<1>::transfer(lst);
}


inline Foam::Xfer<Foam::PackedBoolList> Foam::PackedBoolList::xfer()
{
    return xferMove(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::PackedBoolList&
Foam::PackedBoolList::operator=(const bool val)
{
    PackedList<1>::operator=(val);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator=(const PackedBoolList& lst)
{
    PackedList<1>::operator=(lst);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator=(const PackedList<1>& lst)
{
    PackedList<1>::operator=(lst);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator=(const UList<label>& indices)
{
    clear();
    set(indices);

    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator=(const UIndirectList<label>& indices)
{
    clear();
    set(indices);

    return *this;
}


inline Foam::PackedBoolList
Foam::PackedBoolList::operator~() const
{
    PackedBoolList result(*this);
    result.flip();

    return result;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator&=(const PackedList<1>& lst)
{
    subset(lst);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator&=(const UList<label>& indices)
{
    subset(indices);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator&=(const UIndirectList<label>& indices)
{
    subset(indices);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator|=(const PackedList<1>& lst)
{
    set(lst);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator|=(const UList<label>& indices)
{
    set(indices);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator|=(const UIndirectList<label>& indices)
{
    set(indices);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator+=(const PackedList<1>& lst)
{
    return operator|=(lst);
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator+=(const UList<label>& indices)
{
    return operator|=(indices);
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator+=(const UIndirectList<label>& indices)
{
    return operator|=(indices);
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator-=(const PackedList<1>& lst)
{
    unset(lst);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator-=(const UList<label>& indices)
{
    unset(indices);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator-=(const UIndirectList<label>& indices)
{
    unset(indices);
    return *this;
}


// ************************************************************************* //
