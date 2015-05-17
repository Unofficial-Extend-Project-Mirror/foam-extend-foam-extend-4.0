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
    Foam::atomicWeightTable

Description
    A table of atomic weights for all the elements

SourceFiles
    atomicWeights.C

\*---------------------------------------------------------------------------*/

#ifndef atomicWeights_H
#define atomicWeights_H

#include "scalar.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class atomicWeights Declaration
\*---------------------------------------------------------------------------*/

class atomicWeightTable
:
    public HashTable<scalar>
{

public:

    // Public types and data

        //- Structure to hold the element name and atomic weight pair
        struct atomicWeight
        {
            char name[3];
            scalar weight;
        };

        static const int nElements = 103;

        //- Static table of the weights of all known elements
        static const atomicWeight atomicWeights[nElements];


    // Constructors

        //- Construct from atomicWeights_
        atomicWeightTable();
};


// * * * * * * * * * * * * * * * * Global data  * * * * * * * * * * * * * * //

// Atomic weights table for every element in the periodic table
extern atomicWeightTable atomicWeights;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
