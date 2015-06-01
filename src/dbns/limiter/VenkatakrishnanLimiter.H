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
    VenkatakrishnanLimiter

Description
    Venkatakrishnan differentiable limiter

Author
    Aleksandar Jemcov
    Updated by Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef VenkatakrishnanLimiter_H
#define VenkatakrishnanLimiter_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class VenkatakrishnanLimiter Declaration
\*---------------------------------------------------------------------------*/

class VenkatakrishnanLimiter
{
    // Private data

        //- Limiter value
        scalar k_;

public:

    // Constructor

        //- Construct null
        VenkatakrishnanLimiter()
        :
            k_(5)
        {}

        //- Return limiter
        inline scalar limiter
        (
            const scalar& cellVolume,
            const scalar& deltaOneMax,
            const scalar& deltaOneMin,
            const scalar& deltaTwo
        )
        {
            scalar epsilonSquare = pow3(k_)*cellVolume;

            if (deltaTwo > 0)
            {
                return max
                (
                    min
                    (
                        (
                            (sqr(deltaOneMax) + epsilonSquare)*deltaTwo
                          + 2*sqr(deltaTwo)*deltaOneMax
                        )/
                        stabilise
                        (
                            deltaTwo*
                            (
                                sqr(deltaOneMax)
                              + 2.0*sqr(deltaTwo)
                              + deltaOneMax*deltaTwo
                              + epsilonSquare
                            ),
                            SMALL
                        ),
                        1
                    ),
                    0
                );
            }
            else if (deltaTwo < 0)
            {
                return max
                (
                    min
                    (
                        (
                            (sqr(deltaOneMin) + epsilonSquare)*deltaTwo
                          + 2*sqr(deltaTwo)*deltaOneMin
                        )/
                        stabilise
                        (
                            deltaTwo*
                            (
                                sqr(deltaOneMin)
                              + 2.0*sqr(deltaTwo) + deltaOneMin*deltaTwo
                              + epsilonSquare
                            ),
                            SMALL
                        ),
                        1
                    ),
                    0
                );
            }
            else
            {
                return 1;
            }
        }

        //- Return limiter for a scalar
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
                d & gradPhiP
            );
        }

        //- Return limiter for a vector
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
