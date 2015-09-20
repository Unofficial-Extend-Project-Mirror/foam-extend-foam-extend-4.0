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
    Foam::turbGen

Description
    Generate a turbulent velocity field conforming to a given
    energy spectrum and being divergence free.

SourceFiles
    turbGen.C
    turbGenIO.C

\*---------------------------------------------------------------------------*/

#ifndef turbGen_H
#define turbGen_H

#include "Random.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Kmesh;

/*---------------------------------------------------------------------------*\
                           Class turbGen Declaration
\*---------------------------------------------------------------------------*/

class turbGen
{
    // Private data

        const Kmesh& K;
        const scalar Ea;
        const scalar k0;
        Random RanGen;


public:

    // Constructors

        //- Construct from components
        turbGen(const Kmesh& k, const scalar EA, const scalar K0);


    // Member Functions

        //- Generate and return a velocity field
        vectorField U();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
