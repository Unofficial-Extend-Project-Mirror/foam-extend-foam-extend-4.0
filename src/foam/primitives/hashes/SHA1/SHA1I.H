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

#include "SHA1.H"
#include <cstring>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::SHA1::SHA1()
{
    clear();
}


inline Foam::SHA1::SHA1(const std::string& str)
{
    clear();
    append(str);
}


inline Foam::SHA1::SHA1(const char* str)
{
    clear();
    append(str);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::SHA1& Foam::SHA1::append(const char* data, size_t len)
{
    processBytes(data, len);
    return *this;
}


inline Foam::SHA1& Foam::SHA1::append(const std::string& str)
{
    processBytes(str.data(), str.size());
    return *this;
}


inline Foam::SHA1& Foam::SHA1::append(const char* str)
{
    if (str)
    {
        processBytes(str, strlen(str));
    }
    return *this;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::SHA1::operator==(const SHA1Digest& rhs) const
{
    return this->digest() == rhs;
}


inline bool Foam::SHA1::operator!=(const SHA1Digest& rhs) const
{
    return this->digest() != rhs;
}


inline bool Foam::SHA1::operator==(const SHA1& rhs) const
{
    return digest() == rhs.digest();
}


inline bool Foam::SHA1::operator!=(const SHA1& rhs) const
{
    return digest() != rhs.digest();
}


inline Foam::SHA1::operator
Foam::SHA1Digest () const
{
    return digest();
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Ostream& os, const SHA1& sha)
{
    return os << sha.digest();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
