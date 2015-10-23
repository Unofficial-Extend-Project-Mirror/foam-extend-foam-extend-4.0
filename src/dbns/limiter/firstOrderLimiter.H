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
    firstOrderLimiter

Description
    First order limiter: all second order terms are removed

Author
    Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef firstOrderLimiter_H
#define firstOrderLimiter_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class firstOrderLimiter Declaration
\*---------------------------------------------------------------------------*/

class firstOrderLimiter
{
public:

    // Constructor

        //- Construct null
        firstOrderLimiter()
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
            return 0;
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
            return 0;
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
            return vector::zero;
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
            return symmTensor::zero;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
