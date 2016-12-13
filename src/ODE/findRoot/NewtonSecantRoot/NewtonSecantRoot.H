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
    NewtonSecantRoot

Description
    NewtonSecant root is based on combination of Newton and Secant methods.
    Function and the funciton derivative are provided as a template parameter
    function object, evaluated using operator()(const scalar x)

Author
    Aleksandar Jemcov.  All rights reserved.

SourceFiles
    NewtonSecantRoot.C

\*---------------------------------------------------------------------------*/

#ifndef NewtonSecantRoot_H
#define NewtonSecantRoot_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class NewtonSecantRoot Declaration
\*---------------------------------------------------------------------------*/

template<typename Func, typename Deriv>
class NewtonSecantRoot
{
    // Private data

        //- Reference to a function
        const Func&  f_;

        //- Reference to a derivative
        const Deriv& d_;

        //- Required accuracy
        const scalar eps_;


public:

    // Static data members

        //- Maximum number of iterations
        static const label maxIter;


    // Constructors

        //- Construct given a function
        NewtonSecantRoot(const Func& f, const Deriv& d, const scalar eps);


    // Destructor - default


    // Member Functions

        //- Return root
        scalar root(scalar xN) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NewtonSecantRoot.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
