/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Foam::fieldDictionary

Description
    Read field as dictionary (without mesh).

\*---------------------------------------------------------------------------*/

#ifndef fieldDictionary_H
#define fieldDictionary_H

#include "regIOobject.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class fieldDictionary Declaration
\*---------------------------------------------------------------------------*/

class fieldDictionary
:
    public regIOobject,
    public dictionary
{
    const word type_;

public:

    // Redefine type name to be of the instantiated type
    virtual const word& type() const
    {
        return type_;
    }


    // Constructors

        //- Construct from ioobject and overloaded typename.
        explicit fieldDictionary(const IOobject& io, const word& type)
        :
            regIOobject(io),
            dictionary(readStream(type)),
            type_(type)
        {
            close();
        }



    // Member functions

        bool writeData(Ostream& os) const
        {
            static_cast<const dictionary&>(*this).write(os, false);

            return os.good();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
