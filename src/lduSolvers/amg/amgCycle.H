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
    amgCycle

Description
    Algebraic multigrid fixed cycle class

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    amgCycle.C

\*---------------------------------------------------------------------------*/

#ifndef amgCycle_H
#define amgCycle_H

#include "amgLevel.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class amgCycle Declaration
\*---------------------------------------------------------------------------*/

class amgCycle
{
    // Private data

        //- Pointer to current AMG level
        autoPtr<amgLevel> levelPtr_;

        //- Pointer to coarse AMG cycle
        amgCycle* coarseLevelPtr_;

        //- Number of coarse levels
        label nLevels_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        amgCycle(const amgCycle&);

        //- Disallow default bitwise assignment
        void operator=(const amgCycle&);


public:

    // Static data members

        //- Cycle type
        enum cycleType
        {
            V_CYCLE,
            W_CYCLE,
            F_CYCLE
        };

        //- Static data

        //- Cycle type names
        static const NamedEnum<cycleType, 3> cycleNames_;


    // Constructors

        //- Construct from AMG level
        amgCycle(autoPtr<amgLevel> levelPtr);


    // Destructor

        ~amgCycle();


    // Member Functions

        //- Make coarse levels
        void makeCoarseLevels(const label nMaxLevels);

        //- Return number of levels
        label nLevels() const
        {
            return nLevels_;
        }

        //- Calculate residual
        void residual
        (
            const scalarField& x,
            const scalarField& b,
            const direction cmpt,
            scalarField& res
        ) const
        {
            levelPtr_->residual(x, b, cmpt, res);
        }


        //- Fixed cycle
        void fixedCycle
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt,
            scalarField& xBuffer,
            const cycleType cycle,
            const label nPreSweeps,
            const label nPostSweeps,
            const bool scale
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
