/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Class
    DynList

Description
    A dynamic list is a 1-D vector of objects of type T which resizes
    itself as necessary to accept the new objects.  Internal storage
    is a compact array and the list can be shrunk to compact storage.
    The increase of list size is controlled by three template parameters,
    which allows the list storage to either increase by the given increment
    or the given multiplier and divider (allowing non-integer multiples).

SourceFiles
    DynListI.H
    DynList.C

\*---------------------------------------------------------------------------*/

#ifndef DynList_H
#define DynList_H

#include "UList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class T, label staticSize>
class DynList;

template<class T, label staticSize>
Ostream& operator<<
(
    Ostream&,
    const DynList<T, staticSize>&
);
template<class T, label staticSize>
Istream& operator>>
(
    Istream&,
    DynList<T, staticSize>&
);


/*---------------------------------------------------------------------------*\
                           Class DynList Declaration
\*---------------------------------------------------------------------------*/

template<class T, label staticSize = 16>
class DynList
:
    public UList<T>
{
    // Private data
        //- statically allocated data (used for short lists)
        T staticData_[staticSize];

        //- Number of next free element
        label nextFree_;

    // Private member functions
        //- allocate list size
        inline void allocateSize(const label);

        //- check if index is inside the scope (used for debugging only)
        inline void checkIndex(const label) const;

public:

    // Constructors

        //- Construct null
        inline DynList();

        //- Construct given size
        explicit inline DynList(const label);

        //- Construct from given size and defualt value
        explicit inline DynList(const label, const T&);

        //- Construct from UList. nextFree_ set to size().
        explicit inline DynList(const UList<T>&);

        //- Construct from other ListType
        template<class ListType>
        inline DynList(const ListType&);

        //- Copy constructor
        inline DynList(const DynList<T, staticSize>&);

        //- Construct from Istream. nextFree_ set to size().
        explicit DynList(Istream&);

    // Destructor

        inline ~DynList();


    // Member Functions

        // Access

            //- Size of the active part of the list.
            //- Direct over-ride of list size member function
            inline label size() const;

            //- Number of bytes used by the active part of the list
            //- Direct over-ride of list byteSize member function
            inline label byteSize() const;


        // Edit

            //- Reset size of List.
            void setSize(const label);

            //- Clear the list, i.e. set next free to zero.
            //  Allocated size does not change
            void clear();

            //- Shrink the List<T> to the number of elements used
            void shrink();


    // Member Operators

        //- Append an element at the end of the list
        inline void append(const T& e);

        //- Append an element at the end of the list if it is not yet
        //- present in the list (takes linear time)
        inline void appendIfNotIn(const T& e);

        //- check if the element is in the list (takes linear time)
        inline bool contains(const T& e) const;
        inline label containsAtPosition(const T& e) const;

        //- return a const reference to the last element
        inline const T& lastElement() const;

        //- Return and remove the last element
        inline T removeLastElement();
        inline T removeElement(const label i);

        //- return a refence to the element. Resize the list if necessary
        inline T& newElmt(const label);

        //- Return non-const access to an element,
        //-  resizing the list if necessary
        inline T& operator()(const label);

        //- return access to an element
        inline const T& operator[](const label) const;
        inline T& operator[](const label);

        //- return forward and reverse circular indices
        inline label fcIndex(const label index, const label offset = 1) const;
        inline label rcIndex(const label index, const label offset = 1) const;

        //- return forward and reverse circular elements
        inline const T& fcElement
        (
            const label index,
            const label offset = 1
        ) const;

        inline const T& rcElement
        (
            const label index,
            const label offset = 1
        ) const;

        //- Assignment of all entries to the given value
        inline void operator=(const T&);

        //- Copy of another list
        inline void operator=(const DynList<T, staticSize>&);

        //- Copy of another list type
        template<class ListType>
        inline void operator=(const ListType&);


    // IOstream operators

        // Write DynList to Ostream.
        friend Ostream& operator<< <T, staticSize>
        (
            Ostream&,
            const DynList<T, staticSize>&
        );

        //- Read from Istream, discarding contents of existing DynList.
        friend Istream& operator>> <T, staticSize>
        (
            Istream&,
            DynList<T, staticSize>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DynListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DynList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
