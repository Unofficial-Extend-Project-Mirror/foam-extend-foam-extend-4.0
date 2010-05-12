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
    Foam::UList

Description
    A 1D vector of objects of type \<T\>, where the size of the vector is
    known and used for subscript bounds checking, etc.

    Storage is not allocated during construction or use but is supplied to
    the constructor as an argument.  This type of list is particularly useful
    for lists that refer to parts of existing lists such as SubList.

SourceFiles
    UList.C
    UListI.H
    UListIO.C

\*---------------------------------------------------------------------------*/

#ifndef UList_H
#define UList_H

#include "label.H"
#include "bool.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Pre-declare related List type
template<class T>
class List;

//- Pre-declare related SubList type
template<class T>
class SubList;

// Forward declaration of friend functions and operators

template<class T> class UList;

template<class T> Ostream& operator<<(Ostream&, const UList<T>&);


/*---------------------------------------------------------------------------*\
                           Class UList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UList
{
    // Private data

        //- Vector of values of type T.
        T* __restrict__ v_;

        //- Number of elements in UList.
        label size_;

protected:

    // Protected member functions

        //- Reset the list without changing the base pointer
        //  Used in Sublist slicing.  HJ, 19/Oct/2008
        void reset(T* __restrict__ v, label size);


public:

    // Related types

        //- Declare friendship with the List class
        friend class List<T>;

        //- Declare friendship with the SubList class
        friend class SubList<T>;


    // Constructors

        //- Null constructor.
        inline UList();

        //- Construct from components
        inline UList(T* __restrict__ v, label size);


    // Member functions

        //- Return a null UList
        static const UList<T>& null();


        // Access

            //- Return the forward circular index, i.e. the next index
            //  which returns to the first at the end of the list
            inline label fcIndex(const label i) const;

            //- Return the reverse circular index, i.e. the previous index
            //  which returns to the last at the begining of the list
            inline label rcIndex(const label i) const;

            //- Return the binary size in number of characters of the UList
            //  if the element is a primitive type
            //  i.e. contiguous<T>() == true
            label byteSize() const;


        // Check

            //- Check start is within valid range (0 ... size-1).
            inline void checkStart(const label start) const;

            //- Check size is within valid range (0 ... size).
            inline void checkSize(const label size) const;

            //- Check index i is within valid range (0 ... size-1).
            inline void checkIndex(const label i) const;


        //- Write the UList as a dictionary entry.
        void writeEntry(Ostream& os) const;

        //- Write the UList as a dictionary entry with keyword.
        void writeEntry(const word& keyword, Ostream& os) const;

        //- Assign elements to those from UList.
        void assign(const UList<T>&);


    // Member operators

        //- Return subscript-checked element of UList.
        inline T& operator[](const label);

        //- Return subscript-checked element of constant UList.
        inline const T& operator[](const label) const;

        //- Allow cast to a const List<T>&
        inline operator const Foam::List<T>&() const;

        //- Assignment of all entries to the given value
        void operator=(const T&);


    // STL type definitions

        //- Type of values the UList contains.
        typedef T value_type;

        //- Type that can be used for storing into
        //  UList::value_type objects.
        typedef T& reference;

        //- Type that can be used for storing into
        //  constant UList::value_type objects
        typedef const T& const_reference;

        //- The type that can represent the difference between any two
        //  UList iterator objects.
        typedef label difference_type;

        //- The type that can represent the size of a UList.
        typedef label size_type;


    // STL iterator

        //- Random access iterator for traversing UList.
        typedef T* iterator;

        //- Return an iterator to begin traversing the UList.
        inline iterator begin();

        //- Return an iterator to end traversing the UList.
        inline iterator end();


    // STL const_iterator

        //- Random access iterator for traversing UList.
        typedef const T* const_iterator;

        //- Return a const_iterator to begin traversing the
        //  constant UList.
        inline const_iterator begin() const;

        //- Return a const_iterator to end traversing the
        //  constant UList.
        inline const_iterator end() const;


    // STL reverse_iterator

        //- Reverse iterator for reverse traversal of UList.
        typedef T* reverse_iterator;

        //- Return a reverse_iterator to begin reverse traversing the
        //  UList.
        inline reverse_iterator rbegin();

        //- Return a reverse_iterator to end reverse traversing the
        //  UList.
        inline reverse_iterator rend();


    // STL const_reverse_iterator

        //- Reverse iterator for reverse traversal of constant UList.
        typedef const T* const_reverse_iterator;

        //- Return a const_reverse_iterator to begin reverse traversing the
        //  UList.
        inline const_reverse_iterator rbegin() const;

        //- Return a const_reverse_iterator to end reverse traversing the
        //  UList.
        inline const_reverse_iterator rend() const;


    // STL member functions

        //- Return the number of elements in the UList.
        inline label size() const;

        //- Return size of the largest possible UList.
        inline label max_size() const;

        //- Return true if the UList is empty (i.e., if size() == 0).
        inline bool empty() const;

        //- Swap two ULists of the same type in constant time.
        void swap(UList<T>&);


    // STL member operators

        //- Equality operation on ULists of the same type.
        //  Returns true when the ULists are elementwise equal
        //  (using UList::value_type::operator==).  Takes linear time.
        bool operator==(const UList<T>&) const;

        //- The opposite of the equality operation. Takes linear time.
        bool operator!=(const UList<T>&) const;

        //- Compare two ULists lexicographically. Takes linear time.
        bool operator<(const UList<T>&) const;

        //- Compare two ULists lexicographically. Takes linear time.
        bool operator>(const UList<T>&) const;

        //- Return true if !(a > b). Takes linear time.
        bool operator<=(const UList<T>&) const;

        //- Return true if !(a < b). Takes linear time.
        bool operator>=(const UList<T>&) const;


    // Ostream operator

        // Write UList to Ostream.
        friend Ostream& operator<<
        #ifndef __CINT__
        <T>
        #endif
        (Ostream&, const UList<T>&);
};

// Reverse the first n elements of the list
template<class T>
inline void reverse(UList<T>& ul, const label n);

// Reverse all the elements of the list
template<class T>
inline void reverse(UList<T>& ul);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "UListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/**
 * @def forAll(list, i)
 * Loop across all elements in @a list
 * @par Usage
 * @code
 * forAll(anyList, i)
 * {
 *      statements;
 * }
 * @endcode
 * @sa forAllReverse
*/
/**
 * @def forAllReverse(list, i)
 * Reverse loop across all elements in @a list
 * @par Usage
 * @code
 * forAllReverse(anyList, i)
 * {
 *      statements;
 * }
 * @endcode
 * @sa forAll
*/
#define forAll(list, i) \
    for (Foam::label i=0; i<(list).size(); i++)

#define forAllReverse(list, i) \
    for (Foam::label i=(list).size()-1; i>=0; i--)

/**
 * @def forAllIter(Container, container, iter)
 * Iterate across all elements in the @a container object of type
 * @a Container.
 * @par Usage
 * @code
 * forAll(ContainerType, container, iter)
 * {
 *     statements;
 * }
 * @endcode
 * @sa forAllConstIter
*/
#define forAllIter(Container,container,iter)                                   \
    for                                                                        \
    (                                                                          \
        Container::iterator iter = (container).begin();                        \
        iter != (container).end();                                             \
        ++iter                                                                 \
    )

/**
 * @def forAllConstIter(Container, container, iter)
 * Iterate across all elements in the @a container object of type
 * @a Container with const access.
 * @par Usage
 * @code
 * forAllConstIter(ContainerType, container, iter)
 * {
 *     statements;
 * }
 * @endcode
 * @sa forAllIter
*/
#define forAllConstIter(Container,container,iter)                              \
    for                                                                        \
    (                                                                          \
        Container::const_iterator iter = (container).begin();                  \
        iter != (container).end();                                             \
        ++iter                                                                 \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "UList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
