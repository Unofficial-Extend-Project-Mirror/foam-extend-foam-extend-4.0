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

Description
    lnGrad scheme with limited non-orthogonal correction.

    The limiter is controlled by a coefficient with a value between 0 and 1
    which when 0 switches the correction off and the scheme behaves as
    uncorrectedLnGrad, when set to 1 the full correction is applied and the
    scheme behaves as correctedLnGrad and when set to 0.5 the limiter is
    calculated such that the non-orthogonal contribution does not exceed the
    orthogonal part.

\*---------------------------------------------------------------------------*/

#include "fa.H"
#include "limitedLnGrad.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "correctedLnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * *< * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
limitedLnGrad<Type>::~limitedLnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
limitedLnGrad<Type>::correction
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    GeometricField<Type, faePatchField, edgeMesh> corr =
        correctedLnGrad<Type>(this->mesh()).correction(vf);

    edgeScalarField limiter
    (
        min
        (
            limitCoeff_
           *mag(lnGradScheme<Type>::lnGrad(vf, deltaCoeffs(vf), "orthSnGrad"))
           /(
                (1 - limitCoeff_)*mag(corr)
              + dimensionedScalar("small", corr.dimensions(), SMALL)
            ),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    if (fa::debug)
    {
        Info<< "limitedLnGrad :: limiter min: " << min(limiter.internalField())
            << " max: "<< max(limiter.internalField())
            << " avg: " << average(limiter.internalField()) << endl;
    }

    return limiter*corr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
