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
    MichalakGoochLimiter

Description
    Michalak-Gooch limiter, AIAA-2009-954

Author
    Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef MichalakGoochLimiter_H
#define MichalakGoochLimiter_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class MichalakGoochLimiter Declaration
\*---------------------------------------------------------------------------*/

class MichalakGoochLimiter
{
    // Private data

        //- Coefficient
        const scalar k_;


public:

    // Constructors

        //- Construct null
        MichalakGoochLimiter()
        :
            k_(1.5)
        {}


    // Destructor - default


    // Member functions

        //- Set scalar limiter value
        inline void limiter
        (
            scalar& lim,
            const scalar& cellVolume,
            const scalar& deltaOneMax,
            const scalar& deltaOneMin,
            const scalar& extrapolate
        )
        {
            // If there is no span or extrapolation, do not limit.  Note that
            // deltaOneMax and deltaOneMin are bound by zero and min is
            // negative
            if
            (
                deltaOneMax - deltaOneMin < SMALL
             || mag(extrapolate) < SMALL
            )
            {
                return;
            }

            scalar y;

            if (extrapolate > 0)
            {
                y = deltaOneMax/extrapolate;
            }
            else
            {
                y = deltaOneMin/extrapolate;
            }

            if (y < k_)
            {
                const scalar C2 = (3 - 2*k_)/sqr(k_);
                const scalar C3 = -1/(3*sqr(k_)) - 2/(3*k_)*C2;

                lim = min(lim, y + C2*sqr(y) + C3*pow3(y));
            }
        }

        //- Set Type limiter
        template<class Type>
        inline void limiter
        (
            Type& lim,
            const scalar& cellVolume,
            const Type& deltaOneMax,
            const Type& deltaOneMin,
            const Type& extrapolate
        )
        {
            for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
            {
                limiter
                (
                    lim.component(cmpt),
                    cellVolume,
                    deltaOneMax.component(cmpt),
                    deltaOneMin.component(cmpt),
                    extrapolate.component(cmpt)
                );
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
