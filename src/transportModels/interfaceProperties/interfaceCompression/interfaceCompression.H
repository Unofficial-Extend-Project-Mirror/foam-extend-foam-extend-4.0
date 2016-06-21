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
    Foam::interfaceCompressionLimiter

Description
    Interface compression scheme currently based on the generic limited
    scheme although it does not use the NVD/TVD functions.

SourceFiles
    interfaceCompression.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCompression_H
#define interfaceCompression_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class interfaceCompressionLimiter Declaration
\*---------------------------------------------------------------------------*/

class interfaceCompressionLimiter
{

public:

    interfaceCompressionLimiter(Istream&)
    {}

    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const scalar& phiP,
        const scalar& phiN,
        const vector&,
        const scalar&
    ) const
    {
        // Quadratic compression scheme
        // return min(max(4*min(phiP*(1 - phiP), phiN*(1 - phiN)), 0), 1);

        // Quartic compression scheme
        return
            min
            (
                max
                (
                    1 - max
                    (
                        sqr(1 - 4*phiP*(1 - phiP)),
                        sqr(1 - 4*phiN*(1 - phiN))
                    ),
                    0
                ),
                1
            );
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
