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
    Foam::Callback

Description
    Abstract class to handle automatic call-back registration with the
    CallbackRegistry.  Derive from this class and extend by adding the
    appropriate callback functions.

SourceFiles
    Callback.C

\*---------------------------------------------------------------------------*/

#ifndef Callback_H
#define Callback_H

#include "CallbackRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Callback Declaration
\*---------------------------------------------------------------------------*/

template<class CallbackType>
class Callback
:
    public CallbackRegistry<CallbackType>::link
{
    // Private data

        //- CallbackRegistry reference
        CallbackRegistry<CallbackType>& cbr_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const Callback<CallbackType>&);


public:

    // Constructors

        //- Construct and register with CallbackRegistry
        Callback(CallbackRegistry<CallbackType>&);

        //- Construct as copy
        Callback(const Callback<CallbackType>&);


    // Destructor

        virtual ~Callback();


    // Member functions

        // Registration

            //- Register Callback with CallbackRegistry
            bool checkIn();

            //- Check-out Callback from CallbackRegistry
            bool checkOut();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Callback.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
