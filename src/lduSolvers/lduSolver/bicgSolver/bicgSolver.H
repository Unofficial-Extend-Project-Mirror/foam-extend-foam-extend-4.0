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
    bicgSolver

Description
    Preconditioned Bi-Conjugate Gradient solver with run-time selectable
    preconditioning.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    bicgSolver.C

\*---------------------------------------------------------------------------*/

#ifndef bicgSolver_H
#define bicgSolver_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class bicgSolver Declaration
\*---------------------------------------------------------------------------*/

class bicgSolver
:
    public lduSolver
{
    // Private data

        //- Preconditioner
        autoPtr<lduPreconditioner> preconPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        bicgSolver(const bicgSolver&);

        //- Disallow default bitwise assignment
        void operator=(const bicgSolver&);


public:

    //- Runtime type information
    TypeName("BiCG");


    // Constructors

        //- Construct from matrix components and solver data stream
        bicgSolver
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~bicgSolver()
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
