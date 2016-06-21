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

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::debug::controlSwitches<T>::operator=(const T& rhs)
{
    switchValue_ = rhs;
}


template<class T>
inline const T Foam::debug::controlSwitches<T>::operator&(const T& rhs)
{
    return switchValue_ & rhs;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::operator==
(
    const debug::controlSwitches<T>& lhs,
    const debug::controlSwitches<T>& rhs
)
{
    return lhs() == rhs();
}


template<class T>
inline bool Foam::operator==
(
    const T& lhs,
    const debug::controlSwitches<T>& rhs
)
{
    return lhs == rhs();
}


template<class T>
inline bool Foam::operator==
(
    const debug::controlSwitches<T>& lhs,
    const T& rhs
)
{
    return lhs() == rhs;
}


template<class T>
inline bool Foam::operator!=
(
    const debug::controlSwitches<T>& lhs,
    const debug::controlSwitches<T>& rhs
)
{
    return lhs() != rhs();
}


template<class T>
inline bool Foam::operator!=
(
    const T& lhs,
    const debug::controlSwitches<T>& rhs
)
{
    return lhs != rhs();
}


template<class T>
inline bool Foam::operator!=
(
    const debug::controlSwitches<T>& lhs,
    const T& rhs
)
{
    return lhs() != rhs;
}


template<class T>
inline bool Foam::operator>=
(
    const debug::controlSwitches<T>& lhs,
    const debug::controlSwitches<T>& rhs
)
{
    return lhs() >= rhs();
}


template<class T>
inline bool Foam::operator>=
(
    const T& lhs,
    const debug::controlSwitches<T>& rhs
)
{
    return lhs >= rhs();
}


template<class T>
inline bool Foam::operator>=
(
    const debug::controlSwitches<T>& lhs,
    const T& rhs
)
{
    return lhs() >= rhs;
}


template<class T>
inline bool Foam::operator>
(
    const debug::controlSwitches<T>& lhs,
    const debug::controlSwitches<T>& rhs
)
{
    return lhs() > rhs();
}


template<class T>
inline bool Foam::operator>
(
    const T& lhs,
    const debug::controlSwitches<T>& rhs
)
{
    return lhs > rhs();
}


template<class T>
inline bool Foam::operator>
(
    const debug::controlSwitches<T>& lhs,
    const T& rhs
)
{
    return lhs() > rhs;
}


template<class T>
inline bool Foam::operator<=
(
    const debug::controlSwitches<T>& lhs,
    const debug::controlSwitches<T>& rhs
)
{
    return lhs() <= rhs();
}


template<class T>
inline bool Foam::operator<=
(
    const T& lhs,
    const debug::controlSwitches<T>& rhs
)
{
    return lhs <= rhs();
}


template<class T>
inline bool Foam::operator<=
(
    const debug::controlSwitches<T>& lhs,
    const T& rhs
)
{
    return lhs() <= rhs;
}


template<class T>
inline bool Foam::operator<
(
    const debug::controlSwitches<T>& lhs,
    const debug::controlSwitches<T>& rhs
)
{
    return lhs() < rhs();
}


template<class T>
inline bool Foam::operator<
(
    const T& lhs,
    const debug::controlSwitches<T>& rhs
)
{
    return lhs < rhs();
}


template<class T>
inline bool Foam::operator<
(
    const debug::controlSwitches<T>& lhs,
    const T& rhs
)
{
    return lhs() < rhs;
}


// ************************************************************************* //
