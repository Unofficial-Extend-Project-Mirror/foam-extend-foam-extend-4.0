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

#include <iostream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::string::string()
{}


inline Foam::string::string(const std::string& stdStr)
:
    std::string(stdStr)
{}


// Copy character array
inline Foam::string::string(const char* str)
:
    std::string(str)
{}


// Construct from a given number of characters in a character array
inline Foam::string::string(const char* str, const size_type len)
:
    std::string(str, len)
{}


// Construct from a single character
inline Foam::string::string(const char c)
:
    std::string(1, c)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class String>
inline bool Foam::string::valid(const string& s)
{
    bool iv = false;

    for (const_iterator iter = s.begin(); iter != s.end(); iter++)
    {
        if (!String::valid(*iter))
        {
            iv = true;
            break;
        }
    }

    return !iv;
}


template<class String>
inline bool Foam::string::stripInvalid(string& s)
{
    if (!valid<String>(s))
    {
        register size_type nValid=0;
        iterator iter2 = s.begin();

        for
        (
            const_iterator iter1 = iter2;
            iter1 != const_cast<const string&>(s).end();
            iter1++
        )
        {
            register char c = *iter1;

            if (String::valid(c))
            {
                *iter2 = c;
                ++iter2;
                ++nValid;
            }
        }

        s.resize(nValid);

        return true;
    }

    return false;
}


template<class String>
inline String Foam::string::validate(const string& s)
{
    string ss = s;
    stripInvalid<String>(ss);
    return ss;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::string Foam::string::operator()
(
    const size_type i,
    const size_type n
) const
{
    return substr(i, n);
}


inline Foam::string Foam::string::operator()(const size_type n) const
{
    return substr(0, n);
}


inline Foam::string::hash::hash()
{}


inline Foam::string::size_type Foam::string::hash::operator()
(
    const string& key
) const
{
    register size_type hashVal = 0;

    for (string::const_iterator iter=key.begin(); iter!=key.end(); ++iter)
    {
        hashVal = hashVal<<1 ^ *iter;
    }

    return hashVal;
}


inline Foam::string::size_type Foam::string::hash::operator()
(
    const string& key,
    const size_type tableSize
) const
{
    return ::abs(operator()(key)) % tableSize;
}


// ************************************************************************* //
