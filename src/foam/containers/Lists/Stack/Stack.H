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
#include "dynamicLabelList.H"
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
        dynamicLabelList stack_;

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
