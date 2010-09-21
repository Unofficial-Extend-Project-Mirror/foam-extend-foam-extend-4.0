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
    coarseAmgLevel

Description
    Coarse AMG level stores matrix, x and b locally

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    coarseAmgLevel.C

\*---------------------------------------------------------------------------*/

#ifndef coarseAmgLevel_H
#define coarseAmgLevel_H

#include "amgLevel.H"
#include "amgMatrix.H"
#include "lduPrimitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class coarseAmgLevel Declaration
\*---------------------------------------------------------------------------*/

class coarseAmgLevel
:
    public amgLevel
{
    // Private data

        //- Matrix
        autoPtr<amgMatrix> matrixPtr_;

        //- Solution field
        scalarField x_;

        //- RHS field
        scalarField b_;

        //- Dictionary
        const dictionary& dict_;

        //- AMG policy
        autoPtr<amgPolicy> policyPtr_;

        //- Smoother
        autoPtr<lduSmoother> smootherPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        coarseAmgLevel(const coarseAmgLevel&);

        //- Disallow default bitwise assignment
        void operator=(const coarseAmgLevel&);


public:

        //- Construct from components
        coarseAmgLevel
        (
            autoPtr<amgMatrix> matrixPtr,
            const dictionary& dict,
            const word& policyType,
            const label groupSize,
            const label minCoarseEqns,
            const word& smootherType
        );


    // Destructor

        virtual ~coarseAmgLevel();


    // Member Functions

        //- Return reference to x
        virtual scalarField& x();

        //- Return reference to b
        virtual scalarField& b();

        //- Return reference to dictionary
        const dictionary& dict() const
        {
            return dict_;
        }

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
