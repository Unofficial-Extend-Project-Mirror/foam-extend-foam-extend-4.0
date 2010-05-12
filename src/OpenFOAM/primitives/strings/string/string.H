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

Class
    Foam::string

Description
    A class for handling character strings derived from std::string.

    Strings may contain any characters and therefore are delimited by quotes
    for IO : "any list of characters".

    Used as a base class for word and fileName.

See Also
    Foam::dotFoam() for information about the site/user OpenFOAM
    configuration directory

SourceFiles
    string.C
    stringIO.C

\*---------------------------------------------------------------------------*/

#ifndef string_H
#define string_H

#include "char.H"

#include <string>
#include <cstring>
#include <cstdlib>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class string;
Istream& operator>>(Istream&, string&);
Ostream& operator<<(Ostream&, const string&);
Ostream& operator<<(Ostream&, const std::string&);


/*---------------------------------------------------------------------------*\
                           Class string Declaration
\*---------------------------------------------------------------------------*/

class string
:
    public std::string
{
public:

    // Static data members

        static const char* const typeName;
        static int debug;
        static const string null;


    //- Hashing function class
    class hash
    {

    public:

        inline hash();

        inline size_type operator()(const string& key) const;

        inline size_type operator()
        (
            const string& key,
            const size_type tableSize
        ) const;
    };


    // Constructors

        //- Construct null
        inline string();

        //- Construct from std::string
        inline string(const std::string&);

        //- Construct as copy of character array
        inline string(const char*);

        //- Construct as copy of specified number of characters
        inline string(const char*, const size_type);

        //- Construct from a single character
        inline string(const char);

        //- Construct from Istream
        string(Istream&);


    // Member Functions

        // Access

            //- Count and return the number of a given character in the string
            size_type count(const char) const;

            //- Is this string type valid
            template<class String>
            static inline bool valid(const string&);


        // Edit

            //- Strip invalid characters from the given string
            template<class String>
            static inline bool stripInvalid(string&);

            //- Return a valid String from the given string
            template<class String>
            static inline String validate(const string&);

            //- Replace first occurence of sub-string oldStr with newStr
            //  starting at start
            string& replace
            (
                const string& oldStr,
                const string& newStr,
                size_type start = 0
            );

            //- Replace all occurences of sub-string oldStr with newStr
            //  starting at start
            string& replaceAll
            (
                const string& oldStr,
                const string& newStr,
                size_type start = 0
            );

            //- Expand initial tildes and all occurences of environment variables
            //  Expansion includes:
            //  -# environment variables
            //     - "$VAR", "${VAR}"
            //  -# current directory
            //     - leading "./" : the current directory
            //  -# tilde expansion
            //     - leading "~/" : home directory
            //     - leading "~user" : home directory for specified user
            //     - leading "~OpenFOAM" : site/user OpenFOAM configuration directory
            //
            //  @sa
            //  Foam::dotFoam
            string& expand();

            //- Remove repeated characters returning true if string changed
            bool removeRepeated(const char character);

            //- Return string with repeated characters removed
            string removeRepeated(const char character) const;

            //- Remove trailing character returning true if string changed
            bool removeTrailing(const char character);

            //- Return string with trailing character removed
            string removeTrailing(const char character) const;


    // Member Operators

        //- Return the sub-string from the i-th character for @a n characters
        inline string operator()
        (
            const size_type i,
            const size_type n
        ) const;

        //- Return the sub-string from the first character for @a n characters
        inline string operator()
        (
            const size_type n
        ) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, string&);
        friend Ostream& operator<<(Ostream&, const string&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "stringI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
