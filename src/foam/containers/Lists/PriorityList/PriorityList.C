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
    PriorityList

Description
    Priority list with changing priorities for inserted elements.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "PriorityList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::PriorityList<Type>::sortList()
{
    // Reset indices and sorted indices
    for (label i = 0; i < size_; i++)
    {
        indices_[i] = i;
        sortedIndices_[i] = i;
    }

    for (label i = size_/2 - 1; i >= 0; i--)
    {
        this->bisectionSort(i);
    }

    listSorted_ = true;
}


template<class Type>
void Foam::PriorityList<Type>::bisectionSort(const label startIndex)
{
    label i = startIndex;

    for (;;)
    {
        // Find largest node and left and right children
        label n = i;
        label il = 2*i + 1;
        label ir = il + 1;

        if
        (
            il < size_
         && greater(weights_[indices_[il]], weights_[indices_[n]])
        )
        {
            n = il;
        }

        if
        (
            ir < size_
         && greater(weights_[indices_[ir]], weights_[indices_[n]])
        )
        {
            n = ir;
        }

        // End of bisection
        if (n == i) break;

        // Swap i with largest n
        Foam::Swap(indices_[i], indices_[n]);

        // Swap positions in sorted index list
        sortedIndices_[indices_[i]] = i;
        sortedIndices_[indices_[n]] = n;

        // Update for next position
        i = n;
    }
}


template<class Type>
void Foam::PriorityList<Type>::sortUpwards(const label startIndex)
{
    // Set root
    label i = startIndex;

    label n = (i - 1)/2;

    while (i > 0 && greater(weights_[indices_[i]], weights_[indices_[n]]))
    {
        // Swap node i and n if not in proper order
        Foam::Swap(indices_[i], indices_[n]);

        // Swap positions in sorted index list
        sortedIndices_[indices_[i]] = i;
        sortedIndices_[indices_[n]] = n;

        // Update for next position
        i = n;
        n = (i - 1)/2;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct given capacity
template<class Type>
Foam::PriorityList<Type>::PriorityList(const label capacity)
:
    indices_(capacity),
    weights_(capacity),
    sortedIndices_(capacity),
    size_(capacity),
    listSorted_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::label Foam::PriorityList<Type>::removeHead()
{
    if (!listSorted_)
    {
        this->sortList();
    }

    label maxIndex = indices_[0];

    if (--size_ > 0)
    {
        Foam::Swap(indices_[0], indices_[size_]);
        sortedIndices_[indices_[0]] = 0;
        this->bisectionSort(0);
    }

    return maxIndex;
}


template<class Type>
void Foam::PriorityList<Type>::set
(
    const label i,
    const Type& value
)
{
    weights_[i] = value;

    // List is no longer sorted
    listSorted_ = false;
}


template<class Type>
void Foam::PriorityList<Type>::updateWeight
(
    const label i,
    const Type& newWeight
)
{
    if (!listSorted_)
    {
        this->sortList();
    }

    Type delta = newWeight - weights_[i];

    weights_[i] = newWeight;

    if (greater(delta, pTraits<Type>::zero))
    {
        this->sortUpwards(sortedIndices_[i]);
    }
    else
    {
        this->bisectionSort(sortedIndices_[i]);
    }
}


// ************************************************************************* //
