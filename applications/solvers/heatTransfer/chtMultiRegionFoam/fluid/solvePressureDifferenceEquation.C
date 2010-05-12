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

Description
    Solve pressure difference equation

\*---------------------------------------------------------------------------*/

void solvePressureDifferenceEquation
(
    const label corr,
    const label nCorr,
    const label nNonOrthCorr,
    bool& closedVolume,
    volScalarField& pd,
    const dimensionedScalar& pRef,
    const volScalarField& rho,
    const volScalarField& psi,
    const volScalarField& rUA,
    const volScalarField& gh,
    surfaceScalarField& phi
)
{
    closedVolume = pd.needReference();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvm::ddt(psi, pd)
          + fvc::ddt(psi)*pRef
          + fvc::ddt(psi, rho)*gh
          + fvc::div(phi)
          - fvm::laplacian(rho*rUA, pd)
        );

        //pdEqn.solve();
        if (corr == nCorr-1 && nonOrth == nNonOrthCorr)
        {
            pdEqn.solve(pd.mesh().solver(pd.name() + "Final"));
        }
        else
        {
            pdEqn.solve(pd.mesh().solver(pd.name()));
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi += pdEqn.flux();
        }
    }
}
