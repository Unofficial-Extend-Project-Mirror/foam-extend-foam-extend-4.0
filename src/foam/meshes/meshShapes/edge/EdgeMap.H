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
    Foam::EdgeMap

Description
    Map from edge (expressed as its endpoints) to value

\*---------------------------------------------------------------------------*/

#ifndef EdgeMap_H
#define EdgeMap_H

#include "HashTable.H"
#include "edge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class EdgeMap Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class EdgeMap
:
    public HashTable<T, edge, Hash<edge> >
{

public:

    // Constructors

        //- Construct given initial map size
        EdgeMap(const label size = 128)
        :
            HashTable<T, edge, Hash<edge> >(size)
        {}

        //- Construct from Istream
        EdgeMap(Istream& is)
        :
            HashTable<T, edge, Hash<edge> >(is)
        {}

        //- Construct as copy
        EdgeMap(const EdgeMap<T>& map)
        :
            HashTable<T, edge, Hash<edge> >(map)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
