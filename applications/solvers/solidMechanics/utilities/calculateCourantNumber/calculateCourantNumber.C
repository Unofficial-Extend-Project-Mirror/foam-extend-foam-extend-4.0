/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    calculateCourantNumber

Description
    Simple utility which calculate the Courant number for solid mechanics
    models.

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "constitutiveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
# include "setRootCase.H"
# include "createTime.H"
# include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info<< "\nCalculating Courant number\n" << endl;

  // Calculate Courant number for every face

  // Mechanical properties
  volVectorField U
    (
     IOobject
     (
      "U",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh,
     dimensionedVector("zero", dimLength, vector::zero)
     );
  volSymmTensorField sigma
    (
     IOobject
     (
      "sigma",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh,
     dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
     );
  constitutiveModel rheology(sigma, U);
  volScalarField mu = rheology.mu();
  volScalarField lambda = rheology.lambda();
  volScalarField rho = rheology.rho();
  surfaceScalarField Ef =
      fvc::interpolate(mu*(3*lambda + 2*mu)/(lambda+mu), "E");
  surfaceScalarField nuf = fvc::interpolate(lambda/(2*(lambda+mu)), "nu");
  surfaceScalarField rhof = fvc::interpolate(rho);

  surfaceScalarField waveVelocity =
    Foam::sqrt(Ef*(1 - nuf)/(rhof*(1 + nuf)*(1 - 2*nuf)));

  // Courant number
  scalarField Co =
    waveVelocity.internalField()*runTime.deltaT().value()
    *mesh.surfaceInterpolation::deltaCoeffs().internalField();

  // Calculate required time-step for a Courant number of 1.0
  scalar requiredDeltaT = 1.0 /
    gMax
      (
          mesh.surfaceInterpolation::deltaCoeffs().internalField()
          *waveVelocity.internalField()
          );

  scalar averageCo = gAverage(Co);
  scalar maxCo = gMax(Co);
  scalar averageWaveVel = gAverage(waveVelocity);
  scalar maxWaveVel = gMax(waveVelocity);

  Info<< "\nCourant Number\n\tmean: " << averageCo
      << "\n\tmax: " << maxCo << nl
      << "Wave velocity magnitude\n\tmean " << averageWaveVel
      << "\n\tmax: " << maxWaveVel << nl
      << "Time step required for a maximum Courant number of 1.0 is "
      << requiredDeltaT << endl;

  Info<< "\nEnd\n" << endl;

  return(0);
}


// ************************************************************************* //
