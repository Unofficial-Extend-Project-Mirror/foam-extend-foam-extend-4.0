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
    Foam::postfixedSubRegistry

Description
    Registry of regIOobjects

SourceFiles
    postfixedSubRegistry.C

\*---------------------------------------------------------------------------*/

#ifndef postfixedSubRegistry_H
#define postfixedSubRegistry_H

#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class postfixedSubRegistry Declaration
\*---------------------------------------------------------------------------*/

class postfixedSubRegistry
:
    public objectRegistry
{
    // Private Member Functions

        //- Disallow Copy constructor
        postfixedSubRegistry(const postfixedSubRegistry&);

        //- Disallow default bitwise copy construct and assignment
        void operator=(const postfixedSubRegistry&);


public:

    //- Declare type name for this IOobject
    TypeName("postfixedSubRegistry");


    // Constructors

        //- Construct a sub-registry given an IObject to describe the registry
        //  and an initial estimate for the number of entries
        explicit postfixedSubRegistry
        (
            const IOobject& io,
            const label nIoObjects = 128
        );


    // Destructor

        virtual ~postfixedSubRegistry();


    // Member functions

        // Access

            //- Return mangled fileName
            virtual fileName mangleFileName (const fileName&) const;

        // Edit

            //- Add an regIOobject to registry
            virtual bool checkIn(regIOobject&) const;

            //- Remove an regIOobject from registry
            virtual bool checkOut(regIOobject&) const;

        // Reading

            //- Return true if any of the object's files have been modified
            virtual bool modified() const
            {
                return false;
            }

            //- Read object if modified
            virtual bool readIfModified()
            {
                return true;
            }

        // Writing

            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp
            ) const
            {
                return true;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
