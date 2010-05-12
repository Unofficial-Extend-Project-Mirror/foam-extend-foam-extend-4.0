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
    Foam::PackedList

Description
    List of packed unsigned ints. Gets given the number of bits per item.

SourceFiles
    PackedListI.H
    PackedList.C

\*---------------------------------------------------------------------------*/

#ifndef PackedList_H
#define PackedList_H

#include "labelList.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                        Class PackedListName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(PackedList);


/*---------------------------------------------------------------------------*\
                           Class PackedList Declaration
\*---------------------------------------------------------------------------*/

//- For PackedList
class reference
{
    // private data
    unsigned int& elem_;

    unsigned int mask_;

    label startBit_;

public:

    inline reference(unsigned int& elem, unsigned int mask, label startBit)
    :
        elem_(elem),
        mask_(mask),
        startBit_(startBit)
    {}

    inline void operator=(const unsigned int val)
    {
        unsigned int shiftedMask = mask_ << startBit_;

        unsigned int shiftedVal = val << startBit_;

        elem_ = (elem_ & ~shiftedMask) | shiftedVal;
    }

    inline operator unsigned int () const
    {
        return (elem_ >> startBit_) & mask_;
    }

};


template <int nBits>
class PackedList
:
    private List<unsigned int>
{
    // Private data

        //- Number of nBits entries
        label size_;


    // Private Member Functions

        //- Calculate underlying list size
        inline static label intSize(const label sz);

        //- Calculate index into underlying List.
        inline static label intIndex(const label i);

        //- Check index i is within valid range (0 ... size-1).
        inline void checkIndex(const label i) const;

        //- Check value is representable in nBits
        inline void checkValue(const unsigned int val) const;

public:

    // Constructors

        //- Null constructor
        inline PackedList();

        //- Construct with given size. Note: initializes intList to 0.
        inline PackedList(const label size);

        //- Construct with given size and value for all elements.
        PackedList(const label size, const unsigned int val);

        //- Copy constructor.
        PackedList(const PackedList<nBits>& PList);

        //- Construct from labelList.
        PackedList(const labelList&);

        //- Clone
        inline autoPtr<PackedList<nBits> > clone() const;

    // Member Functions

        // Edit

            //- Reset size of List.
            void setSize(const label);

            //- Clear the list, i.e. set size to zero.
            void clear();

            //- Transfer the contents of the argument List into this List
            //  and annull the argument list.
            void transfer(PackedList<nBits>&);


        // Access

            //- Number of packed elements
            inline label size() const;


            //- Get value at index I
            inline unsigned int get(const label i) const;

            //- Set value at index I. Return true if value changed.
            inline bool set(const label i, const unsigned int val);


    // Member operators

            //- Get value at index i
            inline unsigned int operator[](const label i) const;

            //- Set value at index i. Returns proxy which does actual operation
            inline ::Foam::reference operator[](const label i);

            //- Assignment operator. Takes linear time.
            void operator=(const PackedList<nBits>&);

            //- Assignment of all entries to the given value
            inline void operator=(const unsigned int val);

            //- Return as labelList
            labelList operator()() const;


    // Ostream operator

//        // Write PackedList to Ostream.
//        friend Ostream& operator<< <nBits> (Ostream&, const PackedList<nBits>&);
}
;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "PackedListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PackedList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
