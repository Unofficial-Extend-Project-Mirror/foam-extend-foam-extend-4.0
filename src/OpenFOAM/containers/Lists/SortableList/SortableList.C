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

#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from List
template <class Type>
Foam::SortableList<Type>::SortableList(const List<Type>& values)
:
    List<Type>(values),
    indices_(values.size())
{
    sort();
}


// Construct given size. Sort later on.
template <class Type>
Foam::SortableList<Type>::SortableList(const label size)
:
    List<Type>(size),
    indices_(size)
{}


// Construct given size and initial value. Sort later on.
template <class Type>
Foam::SortableList<Type>::SortableList(const label size, const Type& val)
:
    List<Type>(size, val),
    indices_(size)
{}


// Construct as copy.
template <class Type>
Foam::SortableList<Type>::SortableList(const SortableList<Type>& lst)
:
    List<Type>(lst),
    indices_(lst.indices())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Type>
void Foam::SortableList<Type>::setSize(const label newSize)
{
    List<Type>::setSize(newSize);
    indices_.setSize(newSize);
}


template <class Type>
void Foam::SortableList<Type>::sort()
{
    forAll(indices_, i)
    {
        indices_[i] = i;
    }

    Foam::sort(indices_, less(*this));

    List<Type> tmpValues(this->size());

    forAll(indices_, i)
    {
        tmpValues[i] = this->operator[](indices_[i]);
    }

    List<Type>::transfer(tmpValues);
}



template <class Type>
void Foam::SortableList<Type>::stableSort()
{
    forAll(indices_, i)
    {
        indices_[i] = i;
    }

    Foam::stableSort(indices_, less(*this));

    List<Type> tmpValues(this->size());

    forAll(indices_, i)
    {
        tmpValues[i] = this->operator[](indices_[i]);
    }

    List<Type>::transfer(tmpValues);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template <class Type>
void Foam::SortableList<Type>::operator=(const SortableList<Type>& rhs)
{
    List<Type>::operator=(rhs);
    indices_ = rhs.indices();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
