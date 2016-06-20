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

Class
    Foam::NamedEnum

Description
    Initialise the NamedEnum HashTable from the static list of names.

SourceFiles
    NamedEnum.C

\*---------------------------------------------------------------------------*/

#ifndef NamedEnum_H
#define NamedEnum_H

#include "HashTable.H"
#include "StaticAssert.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Enum, int nEnum>
class NamedEnum;

template<class Enum, int nEnum>
Ostream& operator<<(Ostream&, const NamedEnum<Enum, nEnum>&);


/*---------------------------------------------------------------------------*\
                          Class NamedEnum Declaration
\*---------------------------------------------------------------------------*/

template<class Enum, int nEnum>
class NamedEnum
:
    public HashTable<int>
{
#ifndef SWIG
    //- nEnum must be positive (non-zero)
    StaticAssert(nEnum > 0);
#endif


    // Private Member Functions

        //- Disallow default bitwise copy construct
        NamedEnum(const NamedEnum&);

        //- Disallow default bitwise assignment
        void operator=(const NamedEnum&);


public:

    // Static data members

        //- The set of names corresponding to the enumeration Enum
        static const char* names[nEnum];


    // Constructors

        //- Construct from names
        NamedEnum();


    // Member Functions

        //- Read a word from Istream and return the corresponding
        //  enumeration element
        Enum read(Istream&) const;

        //- Write the name representation of the enumeration to an Ostream
        void write(const Enum e, Ostream&) const;


    // Member Operators

        //- Return the enumeration element corresponding to the given name
        const Enum operator[](const char* name) const
        {
            return Enum(HashTable<int>::operator[](name));
        }

        //- Return the enumeration element corresponding to the given name
        const Enum operator[](const word& name) const
        {
            return Enum(HashTable<int>::operator[](name));
        }

        //- Return the name of the given enumeration element
        const char* operator[](const Enum e) const
        {
            return names[e];
        }

        //- Return the name of the given enumeration element
        //  where e is specified as a int
        const char* operator[](const int e) const
        {
            return names[e];
        }


    // Friend operators

        friend Ostream& operator<< <Enum, nEnum>
        (
            Ostream&,
            const NamedEnum<Enum, nEnum>&
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NamedEnum.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
