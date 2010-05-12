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

#ifndef PackedList_I
#define PackedList_I

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Calculate underlying list size
template<int nBits>
inline label PackedList<nBits>::intSize(const label sz)
{
    const label nElemsPerLabel = sizeof(unsigned int)*8/nBits;

    return (sz+nElemsPerLabel-1)/nElemsPerLabel;
}


// Convert index into index in integer array
template<int nBits>
inline label PackedList<nBits>::intIndex(const label i)
{
    const label nElemsPerLabel = sizeof(unsigned int)*8/nBits;

    // Index in underlying int array
    label elemI = i/nElemsPerLabel;

    return elemI;
}


// Check index i is within valid range (0 ... size-1).
template<int nBits>
inline void PackedList<nBits>::checkIndex(const label i) const
{
    if (!size_)
    {
        FatalErrorIn("PackedList<nBits>::checkIndex(const label)")
            << "attempt to access element from zero sized list"
            << abort(FatalError);
    }
    else if (i<0 || i>=size_)
    {
        FatalErrorIn("PackedList<nBits>::checkIndex(const label)")
            << "index " << i << " out of range 0 ... " << size_-1
            << abort(FatalError);
    }
}


// Check value is representable in nBits
template<int nBits>
inline void PackedList<nBits>::checkValue(const unsigned int val) const
{
    if (val>=(1u << nBits))
    {
        FatalErrorIn("PackedList<T>::set(const unsigned int)")
            << "value " << label(val) << " out of range 0 ... "
            << label((1u << nBits)-1)
            << " representable by " << nBits << " bits"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
template<int nBits>
inline PackedList<nBits>::PackedList()
:
    List<unsigned int>(0),
    size_(0)
{}


// Construct with given size.
template<int nBits>
inline PackedList<nBits>::PackedList(const label size)
:
    List<unsigned int>(intSize(size), 0u),
    size_(size)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int nBits>
inline label PackedList<nBits>::size() const
{
    return size_;
}


// Get value at i
template<int nBits>
inline unsigned int PackedList<nBits>::get(const label i) const
{
#   ifdef DEBUGList
    checkIndex(i);
#   endif

    // Constant: number of elements that fit in an unsigned int
    const label nElemsPerLabel = sizeof(unsigned int)*8/nBits;

    unsigned int mask = ((1u << nBits) - 1);

    label indexInLabel = i % nElemsPerLabel;

    // Starting bit in int.
    label startBit = nBits*indexInLabel;

    return (List<unsigned int>::operator[](intIndex(i)) >> startBit) & mask;
}


template<int nBits>
inline unsigned int PackedList<nBits>::operator[](const label i) const
{
    return get(i);
}


// Set value at i
template<int nBits>
inline bool PackedList<nBits>::set(const label i, const unsigned int val)
{
#   ifdef DEBUGList
    checkIndex(i);
    checkValue(val);
#   endif

    // Constant: number of elements that fit in an unsigned int
    const label nElemsPerLabel = sizeof(unsigned int)*8/nBits;

    unsigned int mask = ((1u << nBits) - 1);

    label indexInLabel = i % nElemsPerLabel;

    // Starting bit in int.
    label startBit = nBits*indexInLabel;


    unsigned int shiftedMask = mask << startBit;

    unsigned int shiftedVal = val << startBit;

    unsigned int& elem = List<unsigned int>::operator[](intIndex(i));

    unsigned int oldElem = elem;

    elem = (elem & ~shiftedMask) | shiftedVal;

    return elem != oldElem;
}


template<int nBits>
inline ::Foam::reference PackedList<nBits>::operator[](const label i)
{
#   ifdef DEBUGList
    checkIndex(i);
#   endif

    // Constant: number of elements that fit in an unsigned int
    const label nElemsPerLabel = sizeof(unsigned int)*8/nBits;

    unsigned int mask = ((1u << nBits) - 1);

    label indexInLabel = i % nElemsPerLabel;

    // Starting bit in int.
    label startBit = nBits*indexInLabel;

    unsigned int& elem = List<unsigned int>::operator[](intIndex(i));

    return ::Foam::reference(elem, mask, startBit);
}


// Set all to val
template<int nBits>
inline void PackedList<nBits>::operator=(const unsigned int val)
{
#   ifdef DEBUGList
    checkValue(val);
#   endif
    List<unsigned int>::operator=(val);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
