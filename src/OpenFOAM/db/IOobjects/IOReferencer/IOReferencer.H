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
    Foam::IOReferencer

Description
    A wrapper class that allows you to retrieve non-IOobjects through the
    object registry.  Performs no input or output, simply for a const reference
    retrieval.

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

        // (Now don't let myObject and myIOobject go out of scope)

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

    const Type * typePtr_;

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
            const Type& reference
        );


    // Destructor

        virtual ~IOReferencer();


    // Member functions

        //- A "do nothing" writeData function, required by regIOobject
        virtual bool writeData(Ostream&) const;
        
        //- Return the pointer as a reference
        const Type& operator()() const;
        
        //- Set the pointer
        void set(const Type& reference);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "IOReferencer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
