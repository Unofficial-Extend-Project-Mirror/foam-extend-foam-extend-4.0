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
    Foam::CallbackRegistry

Description
    Base class with which callbacks are registered.

    Derive from this class and extend by adding the appropriate callback
    functions that loop and call the callback functions for each entry.

SourceFiles
    CallbackRegistry.C

\*---------------------------------------------------------------------------*/

#ifndef CallbackRegistry_H
#define CallbackRegistry_H

#include "UIDLList.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CallbackType>
class Callback;


/*---------------------------------------------------------------------------*\
                       Class CallbackRegistry Declaration
\*---------------------------------------------------------------------------*/

TemplateName(CallbackRegistry);


/*---------------------------------------------------------------------------*\
                       Class CallbackRegistry Declaration
\*---------------------------------------------------------------------------*/

template<class CallbackType>
class CallbackRegistry
:
    public CallbackRegistryName,
    public UIDLList<CallbackType>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        CallbackRegistry(const CallbackRegistry&);

        //- Disallow default bitwise assignment
        void operator=(const CallbackRegistry&);


public:

    // Constructors

        //- Construct null
        CallbackRegistry();


    // Destructor

        virtual ~CallbackRegistry();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CallbackRegistry.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
