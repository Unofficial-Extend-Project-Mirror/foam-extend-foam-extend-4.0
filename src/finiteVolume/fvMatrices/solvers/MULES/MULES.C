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

\*---------------------------------------------------------------------------*/

#include "MULES.H"
#include "upwind.H"
#include "uncorrectedSnGrad.H"
#include "gaussConvectionScheme.H"
#include "gaussLaplacianScheme.H"
#include "uncorrectedSnGrad.H"
#include "surfaceInterpolate.H"
#include "fvcSurfaceIntegrate.H"
#include "slicedSurfaceFields.H"
#include "syncTools.H"

#include "fvCFD.H"

#include "profilingTrigger.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::MULES::explicitSolve
(
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const scalar psiMax,
    const scalar psiMin
)
{
    profilingTrigger trigger("MULES::explicitSolve");
    explicitSolve
    (
        geometricOneField(),
        psi,
        phi,
        phiPsi,
        zeroField(), zeroField(),
        psiMax, psiMin
    );
}


void Foam::MULES::implicitSolve
(
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const scalar psiMax,
    const scalar psiMin
)
{
    profilingTrigger trigger("MULES::implicitSolve");
    implicitSolve
    (
        geometricOneField(),
        psi,
        phi,
        phiPsi,
        zeroField(), zeroField(),
        psiMax, psiMin
    );
}


// ************************************************************************* //
