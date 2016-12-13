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
    Foam::Dictionary

Description
    Gerneral purpose template dictionary class which manages the storage
    associated with it.

    It is derived from DictionaryBase instantiated on a memory managed form
    of intrusive doubly-linked list of \<T\>.

SourceFiles
    DictionaryTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef DictionaryTemplate_H
#define DictionaryTemplate_H

#include "DictionaryBase.H"
#include "IDLList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Dictionary Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Dictionary
:
    public DictionaryBase<IDLList<T>, T>
{

public:

    // Constructors

        //- Null constructor
        Dictionary();

        //- Copy construct
        Dictionary(const Dictionary&);


    // Member functions

        //- Remove an entry specified by keyword and delete the pointer.
        //  Returns true if the keyword was found
        bool erase(const word& keyword);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DictionaryTemplate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
