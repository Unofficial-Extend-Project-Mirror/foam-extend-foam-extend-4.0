/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    PriorityList

Description
    Priority list with changing priorities for inserted elements.
    List is sorted on access to produce element index with highest weight

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

References:
    R. Sedgewick, "Algorithms", Addison-Wesley, Reading, MA, 2nd edition, 1988.

SourceFiles
    PriorityList.C

\*---------------------------------------------------------------------------*/

#ifndef PriorityList_H
#define PriorityList_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class PriorityList Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PriorityList
{
    // Private data

        //- Indices into weights
        labelList indices_;

        //- Weights
        List<Type> weights_;

        //- Sorted indices
        labelList sortedIndices_;

        //- Active size
        label size_;

        //- Is the list sorted?
        mutable bool listSorted_;


    // Private static functions

        //- Greater element comparison
        inline static bool greater(const Type& a, const Type& b)
        {
            return a > b;
        }


    // Private Member Functions

        //- Disallow default bitwise copy construct
        PriorityList(const PriorityList<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const PriorityList<Type>&);


        //- Sort list
        void sortList();

        //- Bisection sort
        void bisectionSort(const label startIndex);

        //- Sort upwards using bisection with root i
        void sortUpwards(const label startIndex);


public:

    // Constructors

        //- Construct given capacity
        PriorityList(const label capacity);


    // Destructor - default


    // Member Functions

        //- Return size
        label size() const
        {
            return size_;
        }

        //- Is the list empty?
        bool empty() const
        {
            return size_ <= 0;
        }

        //- Return indices
        const labelList& indices() const
        {
            return indices_;
        }

        //- Return weights
        const List<Type>& weights() const
        {
            return weights_;
        }


        //- Return index with highest weight and remove
        label removeHead();

        //- Set element with weight
        void set(const label i, const Type& value);

        //- Update weight
        void updateWeight(const label i, const Type& newWeight);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PriorityList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
