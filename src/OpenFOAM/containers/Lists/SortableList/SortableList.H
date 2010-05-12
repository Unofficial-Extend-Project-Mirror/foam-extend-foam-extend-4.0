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
    Foam::SortableList

Description
    A list that is sorted upon construction or when explicitly requested
    with the sort() method.

SourceFiles
    SortableList.C

\*---------------------------------------------------------------------------*/

#ifndef SortableList_H
#define SortableList_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SortableList Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class SortableList
:
    public List<Type>
{
    // Private data

        //- Original indices
        labelList indices_;


public:

    // Public classes

        //- Less function class used by the sort function
        class less
        {
            const UList<Type>& values_;

        public:

            less(const UList<Type>& values)
            :
                values_(values)
            {}

            bool operator()(const label a, const label b)
            {
                return values_[a] < values_[b];
            }
        };


    // Constructors

        //- Construct from List, sorting the elements. Starts with indices set
        //  to index in argument
        explicit SortableList(const List<Type>&);

        //- Construct given size. Sort later on.
        explicit SortableList(const label size);

        //- Construct given size and initial value. Sort later on.
        SortableList(const label size, const Type&);

        //- Construct as copy.
        SortableList(const SortableList<Type>&);


    // Member Functions

        //- Return the list of sorted indices. Updated every sort.
        const labelList& indices() const
        {
            return indices_;
        }

        //- Size the list. If grow can cause undefined indices (until next sort)
        void setSize(const label);

        //- Sort the list (if changed after construction time)
        void sort();

        //- Sort the list (if changed after construction time)
        void stableSort();


    // Member Operators

        void operator=(const SortableList<Type>&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SortableList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
