/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::simpleRegIOobject

Description
    Abstract base class for registered object with I/O. Used in debug symbol
    registration.

SourceFiles
    nullSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef simpleRegIOobject_H
#define simpleRegIOobject_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                         Class simpleRegIOobject Declaration
\*---------------------------------------------------------------------------*/

class simpleRegIOobject
{

public:

    //- Runtime type information
    TypeName("null");


    // Constructors

        //- Construct from objectregistry inserter and name
        simpleRegIOobject
        (
            void (*fn)(const char* name, simpleRegIOobject*),
            const char* name
        )
        {
            (*fn)(name, this);
        }


    //- Destructor
    virtual ~simpleRegIOobject()
    {};


    // Member Functions

        //- Read
        virtual void readData(Istream&) = 0;

        //- Write
        virtual void writeData(Ostream&) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
