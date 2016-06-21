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
    fineAmgLevel

Description
    Fine AMG level refers to matrix, x and b data from the solver or
    preconditioner

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    fineAmgLevel.C

\*---------------------------------------------------------------------------*/

#ifndef fineAmgLevel_H
#define fineAmgLevel_H

#include "amgLevel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fineAmgLevel Declaration
\*---------------------------------------------------------------------------*/

class fineAmgLevel
:
    public amgLevel
{
    // Private data

         //- Reference to matrix
        const lduMatrix& matrix_;

        //- Coupling coefficients, upper
        const FieldField<Field, scalar>& coupleBouCoeffs_;

        //- Coupling coefficients, lower
        const FieldField<Field, scalar>& coupleIntCoeffs_;

        //- Coupling interface fields
        const lduInterfaceFieldPtrsList& interfaceFields_;

        //- Dictionary
        const dictionary& dict_;

        //- AMG policy
        autoPtr<amgPolicy> policyPtr_;

        //- Smoother
        autoPtr<lduSmoother> smootherPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fineAmgLevel(const fineAmgLevel&);

        //- Disallow default bitwise assignment
        void operator=(const fineAmgLevel&);


public:

        //- Construct from components
        fineAmgLevel
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfacesFields,
            const dictionary& dict,
            const word& policyType,
            const label groupSize,
            const label minCoarseEqns,
            const word& smootherType
        );


    // Destructor

        virtual ~fineAmgLevel()
        {}


    // Member Functions

        //- Return reference to dictionary
        const dictionary& dict() const
        {
            return dict_;
        }

        //- Return reference to x
        virtual scalarField& x();

        //- Return reference to b
        virtual scalarField& b();

        //- Calculate residual
        virtual void residual
        (
            const scalarField& x,
            const scalarField& b,
            const direction cmpt,
            scalarField& res
        ) const;

        //- Restrict residual
        virtual void restrictResidual
        (
            const scalarField& x,
            const scalarField& b,
            const direction cmpt,
            scalarField& xBuffer,
            scalarField& coarseRes,
            bool preSweepsDone
        ) const;

        //- Prolongate correction
        virtual void prolongateCorrection
        (
            scalarField& x,
            const scalarField& coarseX
        ) const;

        //- Smooth level
        virtual void smooth
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt,
            const label nSweeps
        ) const;

        //- Solve level
        virtual void solve
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt,
            const scalar tolerance,
            const scalar relTol
        ) const;

        //- Scale x
        virtual void scaleX
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt,
            scalarField& xBuffer
        ) const;

        //- Create next level from current level
        virtual autoPtr<amgLevel> makeNextLevel() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
