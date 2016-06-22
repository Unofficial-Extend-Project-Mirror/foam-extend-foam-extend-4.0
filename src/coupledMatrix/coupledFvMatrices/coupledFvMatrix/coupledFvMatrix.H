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
    coupledFvMatrix

Description
    Coupled Finite-Volume matrix.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    coupledFvMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFvMatrix_H
#define coupledFvMatrix_H

#include "fvMatrices.H"
#include "coupledLduMatrix.H"
#include "className.H"
#include "coupledSolverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class coupledFvMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class coupledFvMatrix
:
    public coupledLduMatrix
{
    // Private member functions

        //- Disallow default bitwise copy construct
        coupledFvMatrix(const coupledFvMatrix<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const coupledFvMatrix<Type>&);


        //- Check size
        void checkSize() const;

        //- Return name for combined variable
        word coupledPsiName() const;


public:

    ClassName("coupledFvMatrix");


    // Constructors

        //- Construct given size
        coupledFvMatrix(const label size)
        :
            coupledLduMatrix(size)
        {}


    // Destructor

        ~coupledFvMatrix()
        {}


    // Member functions

        //- Solve returning the solution statistics.
        //  Solver controls read from dictionary
        coupledSolverPerformance solve(const dictionary&);

        //- Solve returning the solution statistics.
        //  Solver controls read from fvSolution
        coupledSolverPerformance solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "coupledFvMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
