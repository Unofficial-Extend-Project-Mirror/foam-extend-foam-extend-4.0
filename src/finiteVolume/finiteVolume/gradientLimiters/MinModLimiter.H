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
    MinModLimiter

Description
    Non-differentiable cell limiter.  The limiter limits the gradient such
    that the extrapolated face value does not under- or over-shoot the
    neighbouring cell values

Author
    Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef MinModLimiter_H
#define MinModLimiter_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class MinModLimiter Declaration
\*---------------------------------------------------------------------------*/

class MinModLimiter
{
public:

    // Constructors

        //- Construct null
        MinModLimiter()
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
            if (mag(extrapolate) < SMALL)
            {
                lim = 1;
            }
            else
            {
                if (extrapolate > deltaOneMax + VSMALL)
                {
                    lim = min(1, deltaOneMax/extrapolate);
                }
                else if (extrapolate < deltaOneMin - VSMALL)
                {
                    lim = min(1, deltaOneMin/extrapolate);
                }
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
