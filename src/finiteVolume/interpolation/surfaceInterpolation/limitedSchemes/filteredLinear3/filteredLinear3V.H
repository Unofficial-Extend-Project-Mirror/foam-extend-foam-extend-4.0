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
    Foam::filteredLinear3VLimiter

Description
    Class to generate weighting factors for the filteredLinear3V
    differencing scheme.

    The aim is to remove high-frequency modes with "staggering"
    characteristics from vector fields by comparing the face gradient in
    the direction of maximum gradient with both neighbouring cell gradients
    and introduce small amounts of upwind in order to damp these modes.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    filteredLinear3V.C

\*---------------------------------------------------------------------------*/

#ifndef filteredLinear3V_H
#define filteredLinear3V_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class filteredLinear3VLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class filteredLinear3VLimiter
:
    public LimiterFunc
{
    // Private data

        // Scaling corefficient for the gradient ratio,
        // 0 = linear
        // 1 = fully limited
        scalar k_;

public:

    filteredLinear3VLimiter(Istream& is)
    :
        k_(readScalar(is))
    {
        if (k_ < 0 || k_ > 1)
        {
            FatalIOErrorIn("filteredLinear3VLimiter(Istream& is)", is)
                << "coefficient = " << k_
                << " should be >= 0 and <= 1"
                << exit(FatalIOError);
        }
    }

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
        // Difference across face
        vector dfV = phiN - phiP;

        // Scalar difference across the face
        // in the direction in which the difference is largest
        scalar df = dfV & dfV;

        // Twice differences across face-neighbour cells
        // in the direction in which the face-difference is largest
        scalar dP = 2*(dfV & (d & gradcP));
        scalar dN = 2*(dfV & (d & gradcN));

        // Calculate the limiter
        scalar limiter = 1 - k_*(dN - df)*(dP - df)/max(sqr(dN + dP), SMALL);

        // Limit the limiter between linear and upwind
        return max(min(limiter, 1), 0);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
