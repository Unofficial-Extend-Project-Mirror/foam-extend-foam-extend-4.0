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
    Foam::IOReferencer

Description
    A wrapper class that allows you to retrieve non-IOobjects through the
    object registry.  Performs no input or output.  Internally stores pointer
    to object it manages and takes care of its scoping.

    Syntax:
    @verbatim
        // Your non-IOobject, e.g. Type = scalar
        Type myObject;

        // Register it to the object registry
        IOReferencer<scalar> myIOobject
        (
            IOobject
            (
                "lookupName",
                instance,           // can be anything
                [local,]            // (optional) can be anything
                registry,           // can be anything
                IOobject::NO_READ,  // *must* be NO_READ
                IOobject::NO_WRITE  // *must* be NO_WRITE
            ),
            myObject
        );

        // In another part of your solver, e.g. in a custom boundary condition:
        const Type& myConstObject = registry.lookupObject<IOReferencer<Type> >
        (
            "lookupName"
        )();

        // Now you have a local const reference to your non-IOobject
    @endverbatim

SourceFiles
    IOReferencer.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef IOReferencer_H
#define IOReferencer_H

#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class IOReferencer Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IOReferencer
:
    public regIOobject
{
    // Private data

        //- Pointer to object
        Type* typePtr_;


public:

    // Constructors

        //- Construct from IOobject - pointer set to null
        IOReferencer
        (
            const IOobject& io
        );

        //- Construct from IOobject - pointer initialized
        IOReferencer
        (
            const IOobject& io,
            Type* ptr
        );


    // Destructor

        virtual ~IOReferencer();


    // Member functions

        //- A "do nothing" writeData function, required by regIOobject
        virtual bool writeData(Ostream&) const;

        //- Return const reference to object
        const Type& operator()() const;

        //- Return access to object
        Type& operator()();

        //- Set the pointer
        void set(Type* ptr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "IOReferencer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
