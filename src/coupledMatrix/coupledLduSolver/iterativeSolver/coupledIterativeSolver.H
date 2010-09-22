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
    coupledIterativeSolver

Description
    Virtual base class for coupled iterative solvers

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    coupledIterativeSolver.C

\*---------------------------------------------------------------------------*/

#ifndef coupledIterativeSolver_H
#define coupledIterativeSolver_H

#include "coupledLduSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class coupledIterativeSolver Declaration
\*---------------------------------------------------------------------------*/

class coupledIterativeSolver
:
    public coupledLduSolver
{
    // Private data

        //- Control data dictionary
        dictionary dict_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        coupledIterativeSolver(const coupledIterativeSolver&);

        //- Disallow default bitwise assignment
        void operator=(const coupledIterativeSolver&);


protected:

    // Protected data

        //- Solver tolerance
        scalar tolerance_;

        //- Relative tolerance
        scalar relTolerance_;

        //- Minimum number of iterations
        //  (forced irrespective of convergence)
        label minIter_;

        //- Maximum number of iterations
        label maxIter_;


    // Protected Member Functions

        //- Return dictionary
        const dictionary& dict() const;


        //- Is the stop criterion reached
        bool stop(coupledSolverPerformance& solverPerf) const;

        //- Read the control parameters from the dictionary
        virtual void readControls();

        //- Return the matrix norm used to normalise the residual for the
        //  stopping criterion
        scalar normFactor
        (
            const FieldField<Field, scalar>& x,
            const FieldField<Field, scalar>& b,
            const FieldField<Field, scalar>& Ax,
            FieldField<Field, scalar>& tmpField,
            const direction cmpt
        ) const;

        scalar normFactor
        (
            const FieldField<Field, scalar>& x,
            const FieldField<Field, scalar>& b,
            const direction cmpt
        ) const;


public:

    // Constructors

        //- Construct from matrix components and solver data stream
        coupledIterativeSolver
        (
            const word& fieldName,
            const coupledLduMatrix& matrix,
            const PtrList<FieldField<Field, scalar> >& bouCoeffs,
            const PtrList<FieldField<Field, scalar> >& intCoeffs,
            const lduInterfaceFieldPtrsListList& interfaces,
            const dictionary& solverData
        );


    // Destructor

        virtual ~coupledIterativeSolver()
        {}


    // Member functions

        // Access

            scalar tolerance() const
            {
                return tolerance_;
            }

            scalar relTolerance() const
            {
                return relTolerance_;
            }

            label minIter() const
            {
                return minIter_;
            }

            label maxIter() const
            {
                return maxIter_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
