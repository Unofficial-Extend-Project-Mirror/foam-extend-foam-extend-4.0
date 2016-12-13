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

template <typename T>
inline bool operator==(const Foam::safeBool<T>& lhs, bool b)
{
    if (b)
    {
        if (lhs)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        if (lhs)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
}

template <typename T>
inline bool operator==(bool b, const Foam::safeBool<T>& rhs)
{
    if (b)
    {
        if (rhs)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        if (rhs)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
}


template <typename T, typename U>
inline bool operator==
(
    const Foam::safeBool<T>& lhs,
    const Foam::safeBool<U>& rhs
)
{
    lhs.this_type_does_not_support_comparisons();
    return false;
}

template <typename T,typename U>
inline  bool operator!=
(
    const Foam::safeBool<T>& lhs,
    const Foam::safeBool<U>& rhs
)
{
    lhs.this_type_does_not_support_comparisons();
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
