/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2007 H. Jasak All rights reserved
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
    rreAmgSolver

Description
    Reduced Rank Extrapolation Algebraic Multigrid solver with run-time
    selection of policy and cycle

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    rreAmgSolver.C

\*---------------------------------------------------------------------------*/

#ifndef rreAmgSolver_H
#define rreAmgSolver_H

#include "lduMatrix.H"
#include "amgPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class rreAmgSolver Declaration
\*---------------------------------------------------------------------------*/

class rreAmgSolver
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
        rreAmgSolver(const rreAmgSolver&);

        //- Disallow default bitwise assignment
        void operator=(const rreAmgSolver&);


public:

    //- Runtime type information
    TypeName("RREAMG");


    // Constructors

        //- Construct from matrix components and solver data stream
        rreAmgSolver
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~rreAmgSolver()
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
