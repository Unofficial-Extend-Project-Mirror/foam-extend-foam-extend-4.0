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

Application
    hllcALEFlux

Description
    HLLC Arbitrary Lagrangian-Eulerian (ALE) flux class, computes numerical
    fluxes based on HLLC ALE scheme. Uses primitive variables as input and
    gives back conservative numerical fluxes.

    Luo, H.; Baum, Joseph D. and Lohner R. "On the computation of multi-
    material flows using ALE formulation."
    Journal of Computational Physics 194 (2004): 304–328.

Author
    Oliver Borm  All rights reserved.

SourceFiles
    hllcALEFlux.H

\*---------------------------------------------------------------------------*/

#ifndef hllcALEFlux_H
#define hllcALEFlux_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class hllcALEFlux Declaration
\*---------------------------------------------------------------------------*/

class hllcALEFlux
{

public:

    // Constructors

        //- Constructor null
        hllcALEFlux()
        {}


    // Destructor - default


    // Member Functions

        //- Update fluxes based on current left and right states
        void evaluateFlux
        (
            scalar& rhoFlux,
            vector& rhoUFlux,
            scalar& rhoEFlux,
            const scalar& pLeft,
            const scalar& pRight,
            const vector& ULeft,
            const vector& URight,
            const scalar& TLeft,
            const scalar& TRight,
            const scalar& RLeft,
            const scalar& RRight,
            const scalar& CvLeft,
            const scalar& CvRight,
            const vector& Sf,
            const scalar& magSf,
            const vector& dotX
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
