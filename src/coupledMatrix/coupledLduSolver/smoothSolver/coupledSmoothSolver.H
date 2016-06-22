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
    coupledSmoothSolver

Description
    Smoother-solver for coupled matrix

    To improve efficiency, the residual is evaluated after every nSweeps
    smoothing iterations.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    coupledSmoothSolver.C

\*---------------------------------------------------------------------------*/

#ifndef coupledSmoothSolver_H
#define coupledSmoothSolver_H

#include "coupledIterativeSolver.H"
#include "coupledGaussSeidelPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class coupledSmoothSolver Declaration
\*---------------------------------------------------------------------------*/

class coupledSmoothSolver
:
    public coupledIterativeSolver
{
    // Private Data

        //- Number of sweeps before evaluating residual
        label nSweeps_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        coupledSmoothSolver(const coupledSmoothSolver&);

        //- Disallow default bitwise assignment
        void operator=(const coupledSmoothSolver&);


protected:

    // Protected Member Functions

        //- Read the control parameters from the dictionary
        virtual void readControls();

public:

    //- Runtime type information
    TypeName("smoothSolver");


    // Constructors

        //- Construct from matrix
        coupledSmoothSolver
        (
            const word& fieldName,
            const coupledLduMatrix& matrix,
            const PtrList<FieldField<Field, scalar> >& bouCoeffs,
            const PtrList<FieldField<Field, scalar> >& intCoeffs,
            const lduInterfaceFieldPtrsListList& interfaces,
            const dictionary& solverData
        );


    // Destructor

        virtual ~coupledSmoothSolver()
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
