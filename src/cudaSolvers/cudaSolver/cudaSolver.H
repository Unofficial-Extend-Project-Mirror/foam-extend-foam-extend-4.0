/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::cudaSolver

Description
    Base class for GPU based solvers

SourceFiles
    cudaSolver.C

Author
    Dominik Christ, Wikki Ltd.
    Based on Cufflink library by Daniel P. Combest

\*---------------------------------------------------------------------------*/

#ifndef cudaSolver_H
#define cudaSolver_H

#include "cudaTypes.H"
#include "lduMatrix.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class cudaSolver Declaration
\*---------------------------------------------------------------------------*/

class cudaSolver
:
    public lduSolver
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        cudaSolver(const cudaSolver&);

        //- Disallow default bitwise assignment
        void operator=(const cudaSolver&);


public:

    //- Runtime type information
    TypeName("cudaSolver");


    // Constructors

        //- Construct from components
        cudaSolver
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~cudaSolver()
        {}


    // Member Functions
        //- Solve the matrix with this solver
        virtual lduSolverPerformance solve
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt = 0
        ) const = 0;

        //- Create COO matrix from symmetric lduMatrix
        cuspEquationSystem createSymCuspMatrix
        (
            const lduMatrix& matrix,
            const scalarField& x,
            const scalarField& b
        ) const;

        //- Create COO matrix from asymmetric lduMatrix
        cuspEquationSystem createAsymCuspMatrix
        (
            const lduMatrix& matrix,
            const scalarField& x,
            const scalarField& b
        ) const;

        cudaSolverPerformance cudaSolverPerformanceDefault() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
