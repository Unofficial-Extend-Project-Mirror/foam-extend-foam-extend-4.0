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
    Foam::HashSet

Description
    A HashTable with word keys but without contents.

Typedef
    Foam::wordHashSet

Description
    A HashSet with (the default) word keys.

\*---------------------------------------------------------------------------*/

#ifndef HashSet_H
#define HashSet_H

#include "HashTable.H"
#include "empty.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class HashSet Declaration
\*---------------------------------------------------------------------------*/

template<class Key=word, class Hash=string::hash>
class HashSet
:
    public HashTable<empty, Key, Hash>
{

public:

    // Constructors

        //- Construct given initial map size
        explicit HashSet(label size = 100)
        :
            HashTable<empty, Key, Hash>(size)
        {}

        //- Construct from Istream
        explicit HashSet(Istream& is)
        :
            HashTable<empty, Key, Hash>(is)
        {}

        //- Construct from UList of Key
        explicit HashSet(const UList<Key>& ul)
        :
            HashTable<empty, Key, Hash>(2*ul.size())
        {
            forAll (ul, i)
            {
                insert(ul[i]);
            }
        }

        //- Construct as copy
        HashSet(const HashSet<Key, Hash>& hs)
        :
            HashTable<empty, Key, Hash>(hs)
        {}


    // Member Functions

        // Edit

            //- Insert a new hashedEntry
            bool insert(const Key& key)
            {
                return HashTable<empty, Key, Hash>::insert(key, empty());
            }
};


typedef HashSet<> wordHashSet;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
