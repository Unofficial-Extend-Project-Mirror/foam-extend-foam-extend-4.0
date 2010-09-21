/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    deflationSolver

Description
    N-th order deflation solver

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    deflationSolver.C

\*---------------------------------------------------------------------------*/

#ifndef deflationSolver_H
#define deflationSolver_H

#include "lduMatrix.H"
#include "amgPrecon.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class deflationSolver Declaration
\*---------------------------------------------------------------------------*/

class deflationSolver
:
    public lduSolver
{

    typedef FieldField<Field, scalar> scalarFieldField;


    // Private data

        //- Preconditioner
        autoPtr<lduPreconditioner> preconPtr_;

        //- RPM order
        label rpmOrder_;

        //- Maximum number of directions
        label maxDirs_;

        //- Basis tolerance.  Required convergence per iteration
        scalar basisTol_;

        //- Divergence tolerance.  Triggers base increase
        scalar divTol_;

        //- Number of basis steps.  Defines how often subspace
        label nBasisSteps_;

        //- Number of power iterations
        label nPowerIter_;

        //- x buffer
        mutable scalarField xBuffer_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        deflationSolver(const deflationSolver&);

        //- Disallow default bitwise assignment
        void operator=(const deflationSolver&);


        //- Gauss elimination solution on matrix subset
        void GaussSolve
        (
            const label activeSize,
            const scalarSquareMatrix& m,
            scalarField& x,
            const scalarField& b
        ) const;

        //- Calculate residual
        scalar residual
        (
            const scalarField& x,
            const scalarField& b,
            const direction cmpt
        ) const;


public:

    //- Runtime type information
    TypeName("deflation");


    // Constructors

        //- Construct from matrix components and solver data stream
        deflationSolver
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~deflationSolver()
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
