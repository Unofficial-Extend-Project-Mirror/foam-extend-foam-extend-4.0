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
    Solve enthalpy equation

\*---------------------------------------------------------------------------*/

void solveEnthalpyEquation
(
    const volScalarField& rho,
    const volScalarField& DpDt,
    const surfaceScalarField& phi,
    const compressible::RASModel& turb,
    basicThermo& thermo
)
{
    volScalarField& h = thermo.h();

    tmp<fvScalarMatrix> hEqn
    (
        fvm::ddt(rho, h)
      + fvm::div(phi, h)
      - fvm::laplacian(turb.alphaEff(), h)
     ==
        DpDt
    );
    hEqn().relax();
    hEqn().solve();

    thermo.correct();

    Info<< "Min/max T:" << min(thermo.T()) << ' ' << max(thermo.T())
        << endl;
}
