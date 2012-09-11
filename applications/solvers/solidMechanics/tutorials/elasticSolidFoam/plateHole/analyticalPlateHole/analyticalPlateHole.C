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
    Generate analytical solution for a infinite plaste with a circular
    hole.
    Stress field sigma is generated.
    Based on solution outlined in Timoshenko, Theory of Elasticity. 

Author
    plateHoleSolution function by Z. Tukovic
    utility assembled by P. Cardiff

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volFields.H"
#include "fvc.H"
#include "fixedValueFvPatchFields.H"
#include "coordinateSystem.H"

symmTensor plateHoleSolution(const vector& C);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"

#   include "createMesh.H"

  runTime++;

  Info << "Writing analytical solution for an infinite plate with a circular hole,\nwhere"
       << "\n\tradius = 0.5"
       << "\n\tdistant traction = (10,000 0 0 )"
       << nl << endl;
 
  volSymmTensorField sigma
    (
     IOobject
     (
      "analyticalSigma",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
     );

  const volVectorField& C = mesh.C();

  forAll(sigma.internalField(), celli)
    {
      vector curR = vector(C[celli].x(), C[celli].y(), 0);

      sigma.internalField()[celli] = plateHoleSolution(curR);
    }

  forAll(sigma.boundaryField(), patchi)
    {
      forAll(sigma.boundaryField()[patchi], facei)
	{
	  vector curR = vector(C.boundaryField()[patchi][facei].x(), C.boundaryField()[patchi][facei].y(), 0);

	  sigma.boundaryField()[patchi][facei] = plateHoleSolution(curR);
	}
    }

  Info << "Writing analytical sigma tensor" << endl;
  sigma.write();

  Info << nl << "End" << endl;
      
  return 0;
}

// ************************************************************************* //


symmTensor plateHoleSolution(const vector& C)
{
    tensor sigma = tensor::zero;

    scalar T = 10000;
    scalar a = 0.5;

    scalar r = ::sqrt(sqr(C.x()) + sqr(C.y()));
    scalar theta = Foam::atan2(C.y(), C.x());

    coordinateSystem cs("polarCS", C, vector(0, 0, 1), C/mag(C));

    sigma.xx() =
        T*(1 - sqr(a)/sqr(r))/2 
      + T*(1 + 3*pow(a,4)/pow(r,4) - 4*sqr(a)/sqr(r))*::cos(2*theta)/2;

    sigma.xy() =
      - T*(1 - 3*pow(a,4)/pow(r,4) + 2*sqr(a)/sqr(r))*::sin(2*theta)/2;

    sigma.yx() = sigma.xy();

    sigma.yy() =
        T*(1 + sqr(a)/sqr(r))/2 
      - T*(1 + 3*pow(a,4)/pow(r,4))*::cos(2*theta)/2;


    // Transformation to global coordinate system
    sigma = ((cs.R()&sigma)&cs.R().T());

    symmTensor S = symmTensor::zero;

    S.xx() = sigma.xx();
    S.xy() = sigma.xy();
    S.yy() = sigma.yy();
    
    return S;
}
