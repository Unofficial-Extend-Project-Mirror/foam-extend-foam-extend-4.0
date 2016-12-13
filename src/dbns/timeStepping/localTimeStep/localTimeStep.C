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
    Foam::localTimeStep

Author
    Oliver Borm
    Aleksandar Jemcov

\*---------------------------------------------------------------------------*/

#include "localTimeStep.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::localTimeStep::localTimeStep
(
    const volVectorField& U,
    const basicThermo& thermo,
    compressible::turbulenceModel& turbModel
)
:
    mesh_(U.mesh()),
    U_(U),
    thermo_(thermo),
    turbModel_(turbModel),
    CoDeltaT_
    (
        IOobject
        (
            "CoDeltaT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("CoDeltaT", dimTime, 0.1)
    )
{}


void Foam::localTimeStep::update
(
    const scalar maxCo,
    const bool adjustTimeStep
)
{
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    // Compute characteristic length for each cell
    // Calculated from min face delta coefficient.  HJ, 6/Sep/2012
    volScalarField deltaX
    (
        IOobject
        (
            "deltaX",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("great", dimLength, GREAT)
    );

    const surfaceScalarField deltaFace
    (
        1/mesh().surfaceInterpolation::deltaCoeffs()
    );

    // Compute maximum face area for each cell from the internal faces
    forAll (owner, facei)
    {
        deltaX[owner[facei]] =
            min(deltaX[owner[facei]], deltaFace[facei]);

        deltaX[neighbour[facei]] =
            min(deltaX[neighbour[facei]], deltaFace[facei]);
    }

    // Compute maximum face area for each cell from the boundary faces
    forAll (deltaX.boundaryField(), patchi)
    {
        const fvsPatchScalarField& pDeltaFace =
            deltaFace.boundaryField()[patchi];

        const fvPatch& p = pDeltaFace.patch();

        if (p.coupled())
        {
            const unallocLabelList& faceCells = p.patch().faceCells();

            forAll (pDeltaFace, patchFacei)
            {
                deltaX[faceCells[patchFacei]] = min
                (
                    deltaX[faceCells[patchFacei]],
                    pDeltaFace[patchFacei]
                );
            }
        }
    }

    // store kappa as field
    // need a copy here, because the return value is a <tmp>
    const volScalarField kappa = thermo_.Cp()/thermo_.Cv();

    // compute the maximum inviscid deltaT
    volScalarField deltaTInvis
    (
       deltaX/(mag(U_)+sqrt(kappa/thermo_.psi()))
    );

    if (max(turbModel_.muEff()).value() > SMALL)
    {
        // Compute the maximum viscous deltaT
        volScalarField deltaTVis
        (
            sqr(deltaX)*thermo_.rho()/turbModel_.muEff()
        );

        CoDeltaT_ = maxCo*min(deltaTVis, deltaTInvis);
    }
    else
    {
        CoDeltaT_ = maxCo*deltaTInvis;
    }

    if (adjustTimeStep)
    {
        CoDeltaT_ = min(CoDeltaT_);
    }
}


// ************************************************************************* //
