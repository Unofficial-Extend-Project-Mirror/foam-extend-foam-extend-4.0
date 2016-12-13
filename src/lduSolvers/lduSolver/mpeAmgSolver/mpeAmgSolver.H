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
    mpeAmgSolver

Description
    Minimal Polynomial Extrapolation Algebraic Multigrid solver with run-time
    selection of policy and cycle

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    mpeAmgSolver.C

\*---------------------------------------------------------------------------*/

#ifndef mpeAmgSolver_H
#define mpeAmgSolver_H

#include "lduMatrix.H"
#include "amgPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class mpeAmgSolver Declaration
\*---------------------------------------------------------------------------*/

class mpeAmgSolver
:
    public lduSolver
{
private:

    // Private data

        //- AMG preconditioner
        amgPrecon amg_;

        //- Number of smoothing steps before projection
        label nSmoothingSteps_;

        //- Krylov space dimension
        label kDimension_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mpeAmgSolver(const mpeAmgSolver&);

        //- Disallow default bitwise assignment
        void operator=(const mpeAmgSolver&);


public:

    //- Runtime type information
    TypeName("MPEAMG");


    // Constructors

        //- Construct from matrix components and solver data stream
        mpeAmgSolver
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~mpeAmgSolver()
        {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual lduSolverPerformance solve
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
