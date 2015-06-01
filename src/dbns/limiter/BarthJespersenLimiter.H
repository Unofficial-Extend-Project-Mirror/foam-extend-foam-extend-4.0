/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    BarthJespersenLimiter

Description
    Barth-Jespersen limiter

Author
    Aleksandar Jemcov, Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef BarthJespersenLimiter_H
#define BarthJespersenLimiter_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class BarthJespersenLimiter Declaration
\*---------------------------------------------------------------------------*/

class BarthJespersenLimiter
{
public:

    // Constructors

        //- Construct null
        BarthJespersenLimiter()
        {}


    // Desctructor - default


    // Member functions

        //- Return limiter value
        inline scalar limiter
        (
            const scalar& cellVolume,
            const scalar& deltaOneMax,
            const scalar& deltaOneMin,
            const scalar& deltaTwo
        )
        {
            if (mag(deltaTwo) < SMALL)
            {
                return 1;
            }
            else
            {
                scalar deltaTwoStab = stabilise(deltaTwo, VSMALL);

                return min
                (
                    max
                    (
                        max(deltaOneMax/deltaTwoStab, 0),
                        max(deltaOneMin/deltaTwoStab, 0)
                    ),
                    1
                );
            }
        }

        //- Return scalar limiter
        inline scalar limiter
        (
            const scalar& cellVolume,
            const scalar& faceFlux,
            const scalar& deltaOneMax,
            const scalar& deltaOneMin,
            const vector& gradPhiP,
            const vector& gradPhiN,
            const vector& d
        )
        {
            return limiter
            (
                cellVolume,
                deltaOneMax,
                deltaOneMin,
                (d & gradPhiP)
            );
        }

        //- Return vector limiter
        inline vector limiter
        (
            const scalar& cellVolume,
            const scalar& faceFlux,
            const vector& deltaOneMax,
            const vector& deltaOneMin,
            const tensor& gradPhiP,
            const tensor& gradPhiN,
            const vector& d
        )
        {
            vector deltaTwo = d & gradPhiP;

            return vector
            (
                limiter
                (
                    cellVolume,
                    deltaOneMax[vector::X],
                    deltaOneMin[vector::X],
                    deltaTwo[vector::X]
                ),
                limiter
                (
                    cellVolume,
                    deltaOneMax[vector::Y],
                    deltaOneMin[vector::Y],
                    deltaTwo[vector::Y]
                ),
                limiter
                (
                    cellVolume,
                    deltaOneMax[vector::Z],
                    deltaOneMin[vector::Z],
                    deltaTwo[vector::Z]
                )
            );
        }

        //- Return limiter for a symmTensor
        inline symmTensor limiter
        (
            const scalar& cellVolume,
            const scalar& faceFlux,
            const symmTensor& deltaOneMax,
            const symmTensor& deltaOneMin,
            const vector& gradPhiPXX,
            const vector& gradPhiNXX,
            const vector& gradPhiPXY,
            const vector& gradPhiNXY,
            const vector& gradPhiPXZ,
            const vector& gradPhiNXZ,
            const vector& gradPhiPYY,
            const vector& gradPhiNYY,
            const vector& gradPhiPYZ,
            const vector& gradPhiNYZ,
            const vector& gradPhiPZZ,
            const vector& gradPhiNZZ,
            const vector& d
        )
        {
            return symmTensor
            (
                limiter
                (
                    cellVolume,
                    deltaOneMax[0],
                    deltaOneMin[0],
                    d & gradPhiPXX
                ),
                limiter
                (
                    cellVolume,
                    deltaOneMax[1],
                    deltaOneMin[1],
                    d & gradPhiPXY
                ),
                limiter
                (
                    cellVolume,
                    deltaOneMax[2],
                    deltaOneMin[2],
                    d & gradPhiPXZ
                ),
                limiter
                (
                    cellVolume,
                    deltaOneMax[3],
                    deltaOneMin[3],
                    d & gradPhiPYY
                ),
                limiter
                (
                    cellVolume,
                    deltaOneMax[4],
                    deltaOneMin[4],
                    d & gradPhiPYZ
                ),
                limiter
                (
                    cellVolume,
                    deltaOneMax[5],
                    deltaOneMin[5],
                    d & gradPhiPZZ
                )
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
