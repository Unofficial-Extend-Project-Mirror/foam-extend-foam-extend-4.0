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
    Foam::DynamicList

Description
    A 1D vector of objects of type \<T\> which resizes itself as necessary to
    accept the new objects.

    Internal storage is a compact array and the list can be shrunk to compact
    storage. The increase of list size is controlled by three template
    parameters, which allows the list storage to either increase by the given
    increment or the given multiplier and divider (allowing non-integer
    multiples).

SourceFiles
    DynamicListI.H
    DynamicList.C

\*---------------------------------------------------------------------------*/

#ifndef DynamicList_H
#define DynamicList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
class DynamicList;

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Ostream& operator<<
(
    Ostream&,
    const DynamicList<T, SizeInc, SizeMult, SizeDiv>&
);
template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Istream& operator>>
(
    Istream&,
    DynamicList<T, SizeInc, SizeMult, SizeDiv>&
);


/*---------------------------------------------------------------------------*\
                           Class DynamicList Declaration
\*---------------------------------------------------------------------------*/

template<class T, unsigned SizeInc=0, unsigned SizeMult=2, unsigned SizeDiv=1>
class DynamicList
:
    public List<T>
{
    // Private data

        //- Number of next free element
        label nextFree_;


public:

    // Related types

        //- Declare friendship with the List class
        friend class List<T>;

    // Constructors

        //- Construct null
        inline DynamicList();

        //- Construct given size
        explicit inline DynamicList(const label);

        //- Construct from UList. nextFree_ set to size().
        explicit inline DynamicList(const UList<T>&);

        //- Construct from Istream. nextFree_ set to size().
        explicit DynamicList(Istream&);


    // Member Functions

        // Access

            //- Size of the underlying storage.
            inline label allocSize() const;

            //- Size of the active part of the list.
            //  Direct over-ride of list size member function
            inline label size() const;


        // Edit

            //- Reset size of List.
            void setSize(const label);

            //- Reset size of List and value for new elements.
            void setSize(const label, const T&);

            //- Clear the list, i.e. set next free to zero.
            //  Allocated size does not change
            void clear();

            //- Shrink the List<T> to the number of elements used
            inline DynamicList<T, SizeInc, SizeMult, SizeDiv>& shrink();

            //- Transfer the contents of the argument List into this List
            //  and annull the argument list. Is same as List::transfer except
            //  checks that you're not changing the underlying list to something
            //  smaller than nextFree_.
            void transfer(List<T>&);

            //- Transfer the contents of the argument DynamicList into this
            //  DynamicList and annull the argument list.
            void transfer(DynamicList<T, SizeInc, SizeMult, SizeDiv>&);


    // Member Operators

        //- Append an element at the end of the list
        inline void append(const T& e);

        //- Return and remove the top element
        inline T remove();

        //- Return non-const access to an element,
        //  resizing the list if necessary
        inline T& operator()(const label);

        //- Assignment of all entries to the given value
        inline void operator=(const T&);

        //- Assignment to List<T>
        inline void operator=(const List<T>&);


    // IOstream operators

        // Write DynamicList to Ostream.
        friend Ostream& operator<<
        #ifndef __CINT__
        <T, SizeInc, SizeMult, SizeDiv>
        #endif
        (
            Ostream&,
            const DynamicList<T, SizeInc, SizeMult, SizeDiv>&
        );

        //- Read from Istream, discarding contents of existing DynamicList.
        friend Istream& operator>>
        #ifndef __CINT__
        <T, SizeInc, SizeMult, SizeDiv>
        #endif
        (
            Istream&,
            DynamicList<T, SizeInc, SizeMult, SizeDiv>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DynamicListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DynamicList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
