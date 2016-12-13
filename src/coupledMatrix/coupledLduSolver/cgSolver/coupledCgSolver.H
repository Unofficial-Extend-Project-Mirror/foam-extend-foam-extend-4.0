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
    coupledCgSolver

Description
    Preconditioned Conjugate Gradient solver for coupled matrices with
    run-time selectable preconditioning.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    coupledCgSolver.C

\*---------------------------------------------------------------------------*/

#ifndef coupledCgSolver_H
#define coupledCgSolver_H

#include "coupledIterativeSolver.H"
#include "coupledLduPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class coupledCgSolver Declaration
\*---------------------------------------------------------------------------*/

class coupledCgSolver
:
    public coupledIterativeSolver
{
    // Private Data

        // Preconditioner
        autoPtr<coupledLduPrecon> preconPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        coupledCgSolver(const coupledCgSolver&);

        //- Disallow default bitwise assignment
        void operator=(const coupledCgSolver&);


public:

    //- Runtime type information
    TypeName("CG");


    // Constructors

        //- Construct from matrix
        coupledCgSolver
        (
            const word& fieldName,
            const coupledLduMatrix& matrix,
            const PtrList<FieldField<Field, scalar> >& bouCoeffs,
            const PtrList<FieldField<Field, scalar> >& intCoeffs,
            const lduInterfaceFieldPtrsListList& interfaces,
            const dictionary& solverData
        );


    // Destructor

        virtual ~coupledCgSolver()
        {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual coupledSolverPerformance solve
        (
            FieldField<Field, scalar>& x,
            const FieldField<Field, scalar>& b,
            const direction cmpt = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
