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
    amgLevel

Description
    Algebraic multigrid level virtual base class

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    amgLevel.C

\*---------------------------------------------------------------------------*/

#ifndef amgLevel_H
#define amgLevel_H

#include "lduMatrix.H"
#include "amgPolicy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class amgLevel Declaration
\*---------------------------------------------------------------------------*/

class amgLevel
{
public:

    // Destructor

        virtual ~amgLevel()
        {}


    // Member Functions

        //- Return reference to x
        virtual scalarField& x() = 0;

        //- Return reference to b
        virtual scalarField& b() = 0;

        //- Calculate residual
        virtual void residual
        (
            const scalarField& x,
            const scalarField& b,
            const direction cmpt,
            scalarField& res
        ) const = 0;

        //- Restrict residual
        virtual void restrictResidual
        (
            const scalarField& x,
            const scalarField& b,
            const direction cmpt,
            scalarField& xBuffer,
            scalarField& coarseRes,
            bool preSweepsDone
        ) const = 0;

        //- Prolongate correction
        virtual void prolongateCorrection
        (
            scalarField& x,
            const scalarField& coarseX
        ) const = 0;

        //- Smooth level
        virtual void smooth
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt,
            const label nSweeps
        ) const = 0;

        //- Solve level
        virtual void solve
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt,
            const scalar tolerance,
            const scalar relTol
        ) const = 0;

        //- Scale x
        virtual void scaleX
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt,
            scalarField& xBuffer
        ) const = 0;

        //- Create next level from current level
        virtual autoPtr<amgLevel> makeNextLevel() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
