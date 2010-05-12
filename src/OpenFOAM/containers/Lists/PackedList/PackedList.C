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

#include "PackedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct with given size and value for all elements.
template<int nBits>
PackedList<nBits>::PackedList(const label size, const unsigned int val)
:
    List<unsigned int>(intSize(size)),
    size_(size)
{
    for (label i = 0; i < size; i++)
    {
        set(i, val);
    }
}


//- Copy constructor.
template<int nBits>
PackedList<nBits>::PackedList(const PackedList<nBits>& PList)
:
    List<unsigned int>(PList),
    size_(PList.size())
{}


//- Construct from labelList
template<int nBits>
PackedList<nBits>::PackedList(const labelList& lst)
:
    List<unsigned int>(intSize(lst.size()), 0),
    size_(lst.size())
{
    forAll(lst, i)
    {
        set(i, lst[i]);
    }
}


template<int nBits>
autoPtr<PackedList<nBits> > PackedList<nBits>::clone() const
{
    return autoPtr<PackedList<nBits> >(new PackedList<nBits>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <int nBits>
void PackedList<nBits>::setSize(const label size)
{
    List<unsigned int>::setSize(intSize(size));
    size_ = size;
}


template <int nBits>
void PackedList<nBits>::clear()
{
    List<unsigned int>::clear();
    size_ = 0;
}


template <int nBits>
void PackedList<nBits>::transfer(PackedList<nBits>& lst)
{
    size_ = lst.size();
    List<unsigned int>::transfer(lst);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Assignment.
template <int nBits>
void PackedList<nBits>::operator=(const PackedList<nBits>& pl)
{
    setSize(pl.size());
    List<unsigned int>::operator=(pl);
}


template <int nBits>
labelList PackedList<nBits>::operator()() const
{
    labelList elems(size());

    forAll(*this, i)
    {
        elems[i] = get(i);
    }
    return elems;
}


// * * * * * * * * * * * * * * * Ostream Operator *  * * * * * * * * * * * * //

//template <int nBits>
//Ostream& ::Foam::operator<<(Ostream& os, const PackedList<nBits>& PL)
//{
//    os << PL();
//    return os;
//}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
