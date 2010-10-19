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
    Stack

Description
    Member functions of the Stack class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Push items on to the stack
inline void Stack::push(const label index)
{
    stackMutex_.lock();

    if (findIndex(stack_, index) == -1)
    {
        stack_.append(index);
    }

    stackMutex_.unlock();
}


//- Insert item onto stack (no checking)
inline void Stack::insert(const label index)
{
    stack_.append(index);
}


// Pop an item off the stack
inline label Stack::pop()
{
    stackMutex_.lock();

    const label index = stack_.remove();

    stackMutex_.unlock();

    return index;
}


// Remove a specific item off the stack
inline void Stack::remove(const label index)
{
    stackMutex_.lock();

    label loc = findIndex(stack_, index);

    if (loc != -1)
    {
        // Create a new list
        labelList newList(stack_.size() - 1);

        label n = 0;

        // Copy items upto location
        for(label i = 0; i < loc; i++)
        {
            newList[n++] = stack_[i];
        }

        // Copy items from location
        for(label i = (loc + 1); i < stack_.size(); i++)
        {
            newList[n++] = stack_[i];
        }

        // Overwrite
        stack_ = newList;
    }

    stackMutex_.unlock();
}


// Return if the stack is empty or not
inline bool Stack::empty()
{
    return (stack_.size() == 0);
}


//- Return the size of the stack
inline label Stack::size()
{
    return stack_.size();
}


//- Clear out the stack
inline void Stack::clear()
{
    stackMutex_.lock();

    stack_.clear();

    stackMutex_.unlock();
}


//- Print out the stack
inline void Stack::print()
{
    Info << stack_ << endl;
}


//- Access the stack as a list
inline label Stack::operator[]
(
    const label index
)
{
    return stack_[index];
}

} // End namespace Foam

// ************************************************************************* //
