/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Gamma

Description
    Class with operator() which returns the weighting factors for the
    Gamma differencing scheme.  Used in conjunction with the template class
    NVDscheme.

SourceFiles
    GammaMake.C

\*---------------------------------------------------------------------------*/

#ifndef Gamma_H
#define Gamma_H

#include "scalar.H"
#include "vector.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class GammaWeight Declaration
\*---------------------------------------------------------------------------*/

class GammaWeight
{
    scalar k_;

public:

    GammaWeight(Istream& is)
    :
        k_(readScalar(is))
    {
        if (k_ < 0 || k_ > 1)
        {
            FatalIOErrorIn("GammaWeight(Istream& is)", is)
                << "coefficient = " << k_
                << " should be >= 0 and <= 1"
                << exit(FatalIOError);
        }

        // Rescale k_ to be >= 0 and <= 0.5 (TVD conformant)
        // and avoid the /0 when k_ = 0
        k_ = max(k_/2.0, SMALL);
    }


    scalar weight
    (
        scalar cdWeight,
        scalar faceFlux,
        scalar phiP,
        scalar phiN,
        const vector& gradcP,
        const vector& gradcN,
        const vector& d
    ) const
    {
        scalar magd = mag(d);
        vector dHat = d/mag(d);

        scalar gradf = (phiN - phiP)/magd;

        scalar gradcf;
        scalar udWeight;

        if (faceFlux > 0)
        {
            gradcf = dHat & gradcP;
            udWeight = 1;
        }
        else
        {
            gradcf = dHat & gradcN;
            udWeight = 0;
        }

        // Stabilise for division
        gradcf = stabilise(gradcf, SMALL);

        scalar phict = 1 - 0.5*gradf/gradcf;
        scalar limiter = min(max(phict/k_, 0), 1);

        return limiter*cdWeight + (1 - limiter)*udWeight;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
