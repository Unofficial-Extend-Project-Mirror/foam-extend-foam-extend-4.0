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
    Thread-safe stack implementation using an internal DynamicList

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    StackI.H

\*---------------------------------------------------------------------------*/

#ifndef Stack_H
#define Stack_H

#include "label.H"
#include "DynamicList.H"
#include "multiThreader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Stack Declaration
\*---------------------------------------------------------------------------*/

class Stack
{
    // Private data

        //- Internal list
        DynamicList<label> stack_;

        //- Mutex for multi-threading
        Mutex stackMutex_;

public:

    // Constructor
    Stack(){}

    // Member functions for access to the stack

        //- Push items on to the stack
        inline void push(const label index);

        //- Insert item onto stack
        inline void insert(const label index);

        //- Pop an item off the stack
        inline label pop();

        //- Remove a specific index off the stack
        inline void remove(const label index);

        //- Return if a stack is empty or not
        inline bool empty();

        //- Return the size of the stack
        inline label size();

        //- Clear out the stack
        inline void clear();

        //- Print out the stack
        inline void print();

        //- Access the stack as a list
        inline label operator[](const label index);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "StackI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
