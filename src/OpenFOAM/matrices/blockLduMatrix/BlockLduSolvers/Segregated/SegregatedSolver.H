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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    SegregatedSolver

Description
    Segregated solver for block matrices.  Segregation decouples block
    coefficients and calls a scalar solver for each component separately.
    Cross-component coupling is lagged in rhs.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    SegregatedSolver.C

\*---------------------------------------------------------------------------*/

#ifndef SegregatedSolver_H
#define SegregatedSolver_H

#include "blockLduSolvers.H"
#include "BlockIterativeSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class SegregatedSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SegregatedSolver
:
    public BlockLduSolver<Type>
{
    // Private data

        //- Scalar x
        scalarField scalarX_;

        //- Scalar matrix
        blockScalarMatrix scalarMatrix_;

        //- Scalar b
        scalarField scalarB_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        SegregatedSolver(const SegregatedSolver<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const SegregatedSolver<Type>&);


public:

    //- Runtime type information
    TypeName("Segregated");


    // Constructors

        //- Construct from matrix components and solver data stream
        SegregatedSolver
        (
            const word& fieldName,
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


        // Destructor

            virtual ~SegregatedSolver()
            {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual BlockSolverPerformance<Type> solve
        (
            Field<Type>& x,
            const Field<Type>& b
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SegregatedSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
