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
    rusanovFlux

Description
    Rusanov flux, as per Rusanov, V.V.: Calculation of Interaction of
    Non-Steady Shock Waves with Obstacles, J. Comp. Math. Physics,
    1:267:279, 1961.

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak

SourceFiles
    rusanovFlux.H

\*---------------------------------------------------------------------------*/

#ifndef rusanovFlux_H
#define rusanovFlux_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class rusanovFlux Declaration
\*---------------------------------------------------------------------------*/

class rusanovFlux
{

public:

    // Constructors

        //- Construct null
        rusanovFlux()
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
            const scalar& magSf
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
