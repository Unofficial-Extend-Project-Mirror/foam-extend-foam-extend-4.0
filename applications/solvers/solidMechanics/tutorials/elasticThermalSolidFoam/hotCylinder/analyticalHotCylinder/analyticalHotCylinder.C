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
    Generate analytical solution for a thick-walled cylinder with a
    temperature gradient.
    Temperature field T and stress field sigma and generated.
    Based on solution outlined in Timoshenko, Theory of Elasticity. 

Author
    philip.cardiff@ucd.ie

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"

#   include "createMesh.H"

  runTime++;

  Info << "Writing analytical solution for a plain strain cylinder with concentric hole,\nwhere"
       << "\n\tinner radius = 0.5"
       << "\n\touter radius = 0.7"
       << "\n\tinner temperature = 100"
       << "\n\touter temperature = 0"
       << "\n\tinner pressure = 0"
       << "\n\touter pressure = 0"
       << "\n\tE = 200e9"
       << "\n\tu = 0.3"
       << "\n\talpha = 1e-5"
       << nl << endl;
 
  //- inner and outer radii and temperatures
  scalar a = 0.5;
  scalar b = 0.7;
  scalar Ti = 100;
  scalar To = 0;

  //- mechanical and thermal properties
  scalar E = 200e9;
  scalar nu = 0.3;
  scalar alpha = 1e-5;

  //- create T field
  volScalarField T
    (
     IOobject
     (
      "analyticalT",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedScalar("zero", dimTemperature, 0.0)
     );
  
  const volVectorField& C = mesh.C();

  //- radial coordinate
  volScalarField radii =
    C.component(vector::X)*C.component(vector::X) + C.component(vector::Y)*C.component(vector::Y);
  forAll(radii.internalField(), celli)
    {
      radii.internalField()[celli] = ::sqrt(radii.internalField()[celli]);
    }
  forAll(radii.boundaryField(), patchi)
    {
      forAll(radii.boundaryField()[patchi], facei)
	{
	  radii.boundaryField()[patchi][facei] = ::sqrt(radii.boundaryField()[patchi][facei]);
	}
    }

  forAll(T.internalField(), celli)
    {
      const scalar& r = radii[celli];

      T.internalField()[celli] =
	( (Ti-To)/Foam::log(b/a) ) * Foam::log(b/r);
    }

  forAll(T.boundaryField(), patchi)
    {
      forAll(T.boundaryField()[patchi], facei)
	{
	  const scalar& r = radii.boundaryField()[patchi][facei];

	  T.boundaryField()[patchi][facei] =
	    ( (Ti-To)/Foam::log(b/a) ) * Foam::log(b/r);
	}
    }

  //- write temperature file
  Info << "Writing analytical termpature field" << endl;
  T.write();


  //- create sigma field
  volScalarField sigmaR
    (
     IOobject
     (
      "sigmaR",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedScalar("zero", dimForce/dimArea, 0.0)
     );

  forAll(sigmaR.internalField(), celli)
    {
      const scalar& r = radii.internalField()[celli];

      sigmaR.internalField()[celli] =
	( (alpha*E*(Ti-To))/(2*(1-nu)*Foam::log(b/a)) ) *
	(-Foam::log(b/r) -( a*a/(b*b - a*a))*(1 - (b*b)/(r*r))*Foam::log(b/a));
    }

  forAll(sigmaR.boundaryField(), patchi)
    {
      forAll(sigmaR.boundaryField()[patchi], facei)
	{
	  const scalar& r = radii.boundaryField()[patchi][facei];
	  
	  sigmaR.boundaryField()[patchi][facei] =
	    ( (alpha*E*(Ti-To))/(2*(1-nu)*Foam::log(b/a)) ) *
	    ( -Foam::log(b/r) - ( a*a/(b*b - a*a))*(1 - (b*b)/(r*r))*Foam::log(b/a) );
	}
    }

  //- write temperature file
  Info << "\nWriting analytical sigmaR field" << endl;
  sigmaR.write();


  volScalarField sigmaTheta
    (
     IOobject
     (
      "sigmaTheta",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedScalar("zero", dimForce/dimArea, 0.0)
     );

  forAll(sigmaTheta.internalField(), celli)
    {
      const scalar& r = radii.internalField()[celli];

      sigmaTheta.internalField()[celli] =
	( (alpha*E*(Ti-To))/(2*(1-nu)*Foam::log(b/a)) ) *
	(1 -Foam::log(b/r) - ( a*a/(b*b - a*a))*(1 + (b*b)/(r*r))*Foam::log(b/a) );
    }

  forAll(sigmaTheta.boundaryField(), patchi)
    {
      forAll(sigmaTheta.boundaryField()[patchi], facei)
	{
	  const scalar& r = radii.boundaryField()[patchi][facei];
	  
	  sigmaTheta.boundaryField()[patchi][facei] =
	    ( (alpha*E*(Ti-To))/(2*(1-nu)*Foam::log(b/a)) ) *
	    (1 -Foam::log(b/r) - ( a*a/(b*b - a*a))*(1 + (b*b)/(r*r))*Foam::log(b/a) );
	}
    }


  //- write temperature file
  Info << "\nWriting analytical sigmaTheta field" << endl;
  sigmaTheta.write();

  volScalarField sigmaZ
    (
     IOobject
     (
      "sigmaZ",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedScalar("zero", dimForce/dimArea, 0.0)
     );

  forAll(sigmaZ.internalField(), celli)
    {
      //- Timoshenko says this but I am not sure I am not sure the BCs in
      //- the z direction
      // sigmaZ.internalField()[celli] =
      // 	( (alpha*E*(Ti-To))/(2*(1-nu)*Foam::log(b/a)) ) *
      // 	(1 - 2*Foam::log(b/r) - ( 2*a*a/(b*b - a*a))*Foam::log(b/a));

      sigmaZ.internalField()[celli] = 
	0.3*(sigmaR.internalField()[celli] + sigmaTheta.internalField()[celli])
	- E*alpha*(T.internalField()[celli]);
    }

  forAll(sigmaZ.boundaryField(), patchi)
    {
      forAll(sigmaZ.boundaryField()[patchi], facei)
	{
	  //- Timoshenko says this but I am not sure I am not sure the BCs in
	  //- the z direction
	  //sigmaZ.boundaryField()[patchi][facei] =
	  //( (alpha*E*(Ti-To))/(2*(1-nu)*Foam::log(b/a)) ) *
	  //(1 - 2*Foam::log(b/r) - ( 2*a*a/(b*b - a*a))*Foam::log(b/a));

	  //-for general 2-D plain strain problems, the axial stress is given by this:
	  sigmaZ.boundaryField()[patchi][facei] = 
	    nu*(sigmaR.boundaryField()[patchi][facei] + sigmaTheta.boundaryField()[patchi][facei])
	    - E*alpha*(T.boundaryField()[patchi][facei]);
	}
    }


  //- write temperature file
  Info << "\nWriting analytical sigmaZ field" << endl;
  sigmaZ.write();


  //- create analytical sigma tensor

  //- create theta field
  volScalarField theta
    (
     IOobject
     (
      "theta",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh,
     dimensionedScalar("zero", dimless, 0.0)
     );
  forAll(theta.internalField(), celli)
    {
      const scalar& x = mesh.C().internalField()[celli][vector::X];
      const scalar& y = mesh.C().internalField()[celli][vector::Y];

      theta.internalField()[celli] = Foam::atan(y/x);
    }

  forAll(theta.boundaryField(), patchi)
    {
      forAll(theta.boundaryField()[patchi], facei)
	{
	  const scalar& x = mesh.C().boundaryField()[patchi][facei][vector::X];
	  const scalar& y = mesh.C().boundaryField()[patchi][facei][vector::Y];
	  
	  theta.boundaryField()[patchi][facei] = Foam::atan(y/x);
	}
    }

  //- rotation matrix to convert polar stresses to cartesian
  volTensorField rotMat
    (
     IOobject
     (
      "rotMat",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh,
     dimensionedTensor("zero", dimless, tensor::zero)
     );

  forAll(rotMat.internalField(), celli)
    {
      const scalar& t = theta.internalField()[celli];

      rotMat.internalField()[celli] = tensor(::cos(t), ::sin(t), 0,
					     -::sin(t), ::cos(t), 0,
					     0, 0, 1);
    }

  forAll(rotMat.boundaryField(), patchi)
    {
      forAll(rotMat.boundaryField()[patchi], facei)
	{
	  const scalar& t = theta.boundaryField()[patchi][facei];
	  
	  rotMat.boundaryField()[patchi][facei] =  tensor(::cos(t), ::sin(t), 0,
					     -::sin(t), ::cos(t), 0,
					     0, 0, 1);
	}
    }

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

  forAll(sigma.internalField(), celli)
    {
      const scalar& r = sigmaR.internalField()[celli];
      const scalar& t = sigmaTheta.internalField()[celli];
      const scalar& z = sigmaZ.internalField()[celli];

      const tensor& rot = rotMat.internalField()[celli];

      symmTensor sigmaCart(r, 0, 0,
			      t, 0,
			         z);

      sigma.internalField()[celli] =
	    symm(rot.T() & sigmaCart & rot);

      //-for general 2-D plain strain problems, the axial stress is given by this:
      //- (which is not equal to the solution by Timoshenko... hmmmnn)
//       sigma.internalField()[celli][symmTensor::ZZ] = 
// 	0.3*(sigma.internalField()[celli][symmTensor::XX] + sigma.internalField()[celli][symmTensor::YY])
// 	- E*alpha*(T.internalField()[celli]);
    }

  forAll(sigma.boundaryField(), patchi)
    {
      forAll(sigma.boundaryField()[patchi], facei)
	{
	  const scalar& r = sigmaR.boundaryField()[patchi][facei];
	  const scalar& t = sigmaTheta.boundaryField()[patchi][facei];
	  const scalar& z = sigmaZ.boundaryField()[patchi][facei];

	  const tensor& rot = rotMat.boundaryField()[patchi][facei];

	  symmTensor sigmaCart(r, 0, 0,
		                  t, 0,
			             z);
	  sigma.boundaryField()[patchi][facei] =
	    symm(rot.T() & sigmaCart & rot);
	}
    }



  Info << "\nWriting analytical sigma tensor" << endl;
  sigma.write();

  Info << nl << "End" << endl;
      
  return 0;
}


// ************************************************************************* //
