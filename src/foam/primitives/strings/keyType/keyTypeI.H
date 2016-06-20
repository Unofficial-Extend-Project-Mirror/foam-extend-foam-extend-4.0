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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::keyType::keyType()
:
    word(),
    isPattern_(false)
{}


inline Foam::keyType::keyType(const keyType& s)
:
    word(s, false),
    isPattern_(s.isPattern())
{}


inline Foam::keyType::keyType(const word& s)
:
    word(s, false),
    isPattern_(false)
{}


// Construct as copy of string. Expect it to be regular expression
inline Foam::keyType::keyType(const string& s)
:
    word(s, false),
    isPattern_(true)
{}


// Construct as copy of character array
inline Foam::keyType::keyType(const char* s)
:
    word(s, false),
    isPattern_(false)
{}


//- Construct as copy of std::string
inline Foam::keyType::keyType
(
    const std::string& s,
    const bool isPattern
)
:
    word(s, false),
    isPattern_(isPattern)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::keyType::isPattern() const
{
    return isPattern_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline const Foam::keyType& Foam::keyType::operator=(const keyType& s)
{
    // Bypass checking
    string::operator=(s);
    isPattern_ = s.isPattern_;
    return *this;
}


inline const Foam::keyType& Foam::keyType::operator=(const word& s)
{
    word::operator=(s);
    isPattern_ = false;
    return *this;
}


inline const Foam::keyType& Foam::keyType::operator=(const string& s)
{
    // Bypass checking
    string::operator=(s);
    isPattern_ = true;
    return *this;
}


inline const Foam::keyType& Foam::keyType::operator=(const char* s)
{
    // Bypass checking
    string::operator=(s);
    isPattern_ = false;
    return *this;
}


// ************************************************************************* //
