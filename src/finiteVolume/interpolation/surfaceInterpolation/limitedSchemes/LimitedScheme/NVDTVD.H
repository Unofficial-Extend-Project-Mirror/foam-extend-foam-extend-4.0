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
    Foam::NVDTVD

Description
    Foam::NVDTVD

\*---------------------------------------------------------------------------*/

#ifndef NVDTVD_H
#define NVDTVD_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class NVDTVD Declaration
\*---------------------------------------------------------------------------*/

class NVDTVD
{

public:

    typedef scalar phiType;
    typedef vector gradPhiType;


    // Null Constructor

        NVDTVD()
        {}


    // Member Functions

        scalar phict
        (
            const scalar faceFlux,
            const scalar phiP,
            const scalar phiN,
            const vector& gradcP,
            const vector& gradcN,
            const vector& d
        ) const
        {
            scalar gradf = phiN - phiP;

            scalar gradcf;

            if (faceFlux > 0)
            {
                gradcf = d & gradcP;
            }
            else
            {
                gradcf = d & gradcN;
            }

            // Stabilise for division
            gradcf = stabilise(gradcf, SMALL);

            return 1 - 0.5*gradf/gradcf;
        }


        scalar r
        (
            const scalar faceFlux,
            const scalar phiP,
            const scalar phiN,
            const vector& gradcP,
            const vector& gradcN,
            const vector& d
        ) const
        {
            scalar gradf = phiN - phiP;

            scalar gradcf;

            if (faceFlux > 0)
            {
                gradcf = d & gradcP;
            }
            else
            {
                gradcf = d & gradcN;
            }

            // Stabilise for division
            // Changed to SMALL to prevent FPE.  OB, 14/Jan/2011
            gradf = stabilise(gradf, SMALL);

            // New formulation.  Oliver Borm and Aleks Jemcov
            // HJ, 13/Jan/2011
            return max(2*(gradcf/gradf) - 1, 0);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
