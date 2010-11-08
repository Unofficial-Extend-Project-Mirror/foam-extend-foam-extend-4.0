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
    BlockGMRESSolver

Description
    Preconditioned Generalised Minimal Residual solver with
    run-time selectable preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    BlockGMRESSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockGMRESSolver_H
#define BlockGMRESSolver_H

#include "blockLduSolvers.H"
#include "BlockIterativeSolver.H"
#include "blockLduPrecons.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class BlockGMRESSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockGMRESSolver
:
    public BlockIterativeSolver<Type>
{
    // Private data

        //- Preconditioner
        autoPtr<BlockLduPrecon<Type> > preconPtr_;

        //- Krylov space dimension
        label nDirs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockGMRESSolver(const BlockGMRESSolver<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const BlockGMRESSolver<Type>&);

        //- Givens rotation
        void givensRotation
        (
            const scalar& H,
            const scalar& beta,
            scalar& c,
            scalar& s
        ) const;


public:

    //- Runtime type information
    TypeName("GMRES");


    // Constructors

        //- Construct from matrix components and solver data stream
        BlockGMRESSolver
        (
            const word& fieldName,
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


        // Destructor

            virtual ~BlockGMRESSolver()
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
#   include "BlockGMRESSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
