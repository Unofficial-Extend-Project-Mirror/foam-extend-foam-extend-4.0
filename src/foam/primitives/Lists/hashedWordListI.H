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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::hashedWordList::contains(const word& name) const
{
    return indices_.found(name);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::hashedWordList::operator=(const UList<word>& lst)
{
    List<word>::operator=(lst);
    rehash();
}


inline void Foam::hashedWordList::operator=(const hashedWordList& lst)
{
    operator=(static_cast<const UList<word>&>(lst));
}


inline const Foam::word& Foam::hashedWordList::operator[]
(
    const label index
) const
{
    return List<word>::operator[](index);
}


// could return -1 instead of bombing out
inline Foam::label Foam::hashedWordList::operator[](const word& name) const
{
    return indices_[name];
}


// ************************************************************************* //
