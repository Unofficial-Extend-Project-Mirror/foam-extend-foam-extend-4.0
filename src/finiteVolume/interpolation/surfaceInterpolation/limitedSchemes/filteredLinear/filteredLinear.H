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
    Foam::filteredLinearLimiter

Description
    Class to generate weighting factors for the filtered-linear
    differencing scheme.

    The aim is to remove high-frequency modes with "staggering"
    characteristics by comparing the face gradient relative to the
    background distribution represented by the neighbouring cell gradients
    with those gradients and introduce small amounts of upwind in order to
    damp these modes.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    filteredLinear.C

\*---------------------------------------------------------------------------*/

#ifndef filteredLinear_H
#define filteredLinear_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class filteredLinearWeight Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class filteredLinearLimiter
:
    public LimiterFunc
{

public:

    filteredLinearLimiter(Istream&)
    {}

    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const typename LimiterFunc::phiType& phiP,
        const typename LimiterFunc::phiType& phiN,
        const typename LimiterFunc::gradPhiType& gradcP,
        const typename LimiterFunc::gradPhiType& gradcN,
        const vector& d
    ) const
    {
        scalar df = phiN - phiP;

        scalar dcP = d & gradcP;
        scalar dcN = d & gradcN;

        scalar limiter =
            2
          - 0.5*min(mag(df - dcP), mag(df - dcN))
           /(max(mag(dcP), mag(dcN)) + SMALL);

        // Limit the limiter between linear and 20% upwind
        return max(min(limiter, 1), 0.8);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
