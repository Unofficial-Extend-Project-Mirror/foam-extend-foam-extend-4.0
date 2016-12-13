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
    Foam::Random

Description
    Simple random number generator.

SourceFiles
    Random.C

\*---------------------------------------------------------------------------*/

#ifndef Random_H
#define Random_H

#include "vector.H"
#include "tensor.H"
#include "symmTensor4thOrder.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Random Declaration
\*---------------------------------------------------------------------------*/

class Random
{
    // Private data

        //- Seed for random numbers
        label Seed;


public:

    // Constructors

        //- Construct given seed
        Random(const label&);


    // Member functions

        //- Return random bit
        int bit();

        //- scalar [0..1] (so including 0,1)
        scalar scalar01();

        //- vector with every component scalar01
        vector vector01();

        //- sphericalTensor with every component scalar01
        sphericalTensor sphericalTensor01();

        //- symmTensor with every component scalar01
        symmTensor symmTensor01();

        //- symmTensor4thOrder with every component scalar01
        symmTensor4thOrder symmTensor4thOrder01();

        //- diagTensor with every component scalar01
         diagTensor diagTensor01();

        //- tensor with every component scalar01
        tensor tensor01();

        //- label [lower..upper]
        label integer(const label lower, const label upper);

        vector position(const vector&, const vector&);

        void randomise(scalar&);
        void randomise(vector&);
        void randomise(sphericalTensor&);
        void randomise(symmTensor&);
        void randomise(symmTensor4thOrder&);
        void randomise(diagTensor&);
        void randomise(tensor&);

        //- Return a normal Gaussian random number
        //  with zero mean and unit variance N(0, 1)
        scalar GaussNormal();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
