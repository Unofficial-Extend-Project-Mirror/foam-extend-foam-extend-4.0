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
    Foam::StaticHashTable

Description
    STL conforming hash table.

    Uses straight lists as underlying type.
    Is slower to insert than the standard HashTable, but should be more
    memory efficient and faster to access.
    Explicitly does not have default size.

SourceFiles
    StaticHashTableI.H
    StaticHashTable.C
    StaticHashTableIO.C

\*---------------------------------------------------------------------------*/

#ifndef StaticHashTable_H
#define StaticHashTable_H

#include "label.H"
#include "word.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T>
class List;

template<class T, class Key, class Hash> class StaticHashTable;

template<class T, class Key, class Hash> Istream& operator>>
(
    Istream&,
    StaticHashTable<T, Key, Hash>&
);

template<class T, class Key, class Hash> Ostream& operator<<
(
    Ostream&,
    const StaticHashTable<T, Key, Hash>&
);


/*---------------------------------------------------------------------------*\
                        Class StaticHashTableName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(StaticHashTable);


/*---------------------------------------------------------------------------*\
                          Class StaticHashTable Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Key=word, class Hash=string::hash>
class StaticHashTable
:
    public StaticHashTableName
{
    // Private data type for table entries

        //- The lookup keys, ordered per hash value
        List<List<Key> > keys_;

        //- For each key the corresponding object.
        List<List<T> > objects_;

        //- The current number of elements in table
        label nElmts_;

public:


    // Forward declaration of STL iterators

        template<class TRef, class TableRef>
        class Iterator;

        typedef Iterator
        <
            T&,
            StaticHashTable<T, Key, Hash>&
        > iterator;

        typedef Iterator
        <
            const T&,
            const StaticHashTable<T, Key, Hash>&
        > const_iterator;


    // Declare friendship with the iterators

        friend class Iterator
        <
            T&,
            StaticHashTable<T, Key, Hash>&
        >;

        friend class Iterator
        <
            const T&,
            const StaticHashTable<T, Key, Hash>&
        >;


    // Constructors

        //- Construct given initial table size
        StaticHashTable(const label size = 100);

        //- Construct from Istream
        StaticHashTable(Istream&, const label size = 100);

        //- Construct as copy
        StaticHashTable(const StaticHashTable<T, Key, Hash>&);


    // Destructor

        ~StaticHashTable();


    // Member Functions

        // Access

            //- Return number of elements in table.
            inline label size() const;

            //- Return true if hashedEntry is found in table
            bool found(const Key& key) const;

            //- Find and return an iterator set at the hashedEntry
            //  If not found iterator = end()
            iterator find(const Key& key);

            //- Find and return an const_iterator set at the hashedEntry
            //  If not found iterator = end()
            const_iterator find(const Key& key) const;

            //- Return the table of contents
            List<Key> toc() const;


        // Edit

            //- Insert a new hashedEntry
            bool insert(const Key& key, const T& newElmt);

            //- Erase an hashedEntry specified by given iterator
            bool erase(const iterator& it);

            //- Erase an hashedEntry specified by given key if in table
            bool erase(const Key& key);

            //- Resize the hash table for efficiency
            void resize(const label newSize);

            //- Clear all entries from table
            void clear();

            //- Transfer the contents of the argument table into this table
            //  and annull the argument table.
            void transfer(StaticHashTable<T, Key, Hash>&);


    // Member Operators

        //- Find and return an hashedEntry
        inline T& operator[](const Key& key);

        //- Find and return an hashedEntry
        inline const T& operator[](const Key& key) const;

        //- Assignment
        void operator=(const StaticHashTable<T, Key, Hash>&);


    // STL type definitions

        //- Type of values the StaticHashTable contains.
        typedef T value_type;

        //- Type that can be used for storing into StaticHashTable::value_type
        //  objects.  This type is usually List::value_type&.
        typedef T& reference;

        //- Type that can be used for storing into constant
        //  StaticHashTable::value_type objects.  This type is usually const
        //  StaticHashTable::value_type&.
        typedef const T& const_reference;

        //- The type that can represent the size of a StaticHashTable.
        typedef label size_type;


    // STL iterator

        //- An STL iterator
        template<class TRef, class TableRef>
        class Iterator
        {
            friend class StaticHashTable;

#           ifndef __INTEL_COMPILER
            template<class TRef2, class TableRef2>
            friend class Iterator;
#           endif

            // Private data

                //- Reference to the StaticHashTable this is an iterator for
                TableRef curStaticHashTable_;

                //- Current hash index
                label hashIndex_;

                //- Index of current element at hashIndex
                label elementIndex_;

        public:

            // Constructors

                //- Construct from hash table, element and hash index
                inline Iterator
                (
                    TableRef curStaticHashTable,
                    label hashIndex_,
                    label elementIndex_
                );

                //- Construct from the non-const iterator
                inline Iterator(const iterator&);


            // Member operators

                inline void operator=(const iterator& iter);

                inline bool operator==(const iterator& iter) const;
                inline bool operator==(const const_iterator& iter) const;

                inline bool operator!=(const iterator& iter) const;
                inline bool operator!=(const const_iterator& iter) const;

                inline TRef operator*();
                inline TRef operator()();

                inline Iterator& operator++();
                inline Iterator operator++(int);

                inline const Key& key();
        };


        //- iterator set to the begining of the StaticHashTable
        inline iterator begin();

        //- iterator set to beyond the end of the StaticHashTable
        inline const iterator& end();

        //- const_iterator set to the begining of the StaticHashTable
        inline const_iterator begin() const;

        //- const_iterator set to beyond the end of the StaticHashTable
        inline const const_iterator& end() const;


    // IOstream Operator

        friend Istream& operator>> <T, Key, Hash>
        (
            Istream&,
            StaticHashTable<T, Key, Hash>&
        );

        friend Ostream& operator<< <T, Key, Hash>
        (
            Ostream&,
            const StaticHashTable<T, Key, Hash>&
        );


private:

        //- iterator returned by end()
        iterator endIter_;

        //- const_iterator returned by end()
        const_iterator endConstIter_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "StaticHashTableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef NoStaticHashTableC
#ifdef NoRepository
#   include "StaticHashTable.C"
#endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
