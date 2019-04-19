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
    Foam::simpleObjectRegistry

Description
    Object registry for simpleRegIOobject. Maintains ordering.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef simpleObjectRegistry_H
#define simpleObjectRegistry_H

#include "Dictionary.H"
#include "simpleRegIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class simpleObjectRegistryEntry Declaration
\*---------------------------------------------------------------------------*/

class simpleObjectRegistryEntry
:
    public Dictionary<simpleObjectRegistryEntry>::link,
    public List<simpleRegIOobject*>
{
public:

    simpleObjectRegistryEntry(const List<simpleRegIOobject*>& data)
    :
        List<simpleRegIOobject*>(data)
    {}
};


/*---------------------------------------------------------------------------*\
                      Class simpleObjectRegistry Declaration
\*---------------------------------------------------------------------------*/

class simpleObjectRegistry
:
    public Dictionary<simpleObjectRegistryEntry>
{
public:

    // Constructors

        //- Construct given initial table size
        simpleObjectRegistry(const label nIoObjects = 128)
        :
            Dictionary<simpleObjectRegistryEntry>(nIoObjects)
        {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
