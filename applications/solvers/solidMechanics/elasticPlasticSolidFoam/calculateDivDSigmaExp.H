if(divDSigmaExpMethod == "standard")
  {
    divDSigmaExp = fvc::div
      (
       mu*gradDU.T() + lambda*(I*tr(gradDU)) - (mu + lambda)*gradDU,
       "div(sigma)"
       );
  }
 else if(divDSigmaExpMethod == "surface")
   { 
     divDSigmaExp = fvc::div
       (
	muf*(mesh.Sf() & fvc::interpolate(gradDU.T()))
	+ lambdaf*(mesh.Sf() & I*fvc::interpolate(tr(gradDU)))
	- (muf + lambdaf)*(mesh.Sf() & fvc::interpolate(gradDU))
	);
   }
 else if(divDSigmaExpMethod == "decompose")
   {
     surfaceTensorField shearGradU =
       ((I - n*n)&fvc::interpolate(gradDU));
     
     divDSigmaExp = fvc::div
       (
	mesh.magSf()
	*(
	  - (muf + lambdaf)*(fvc::snGrad(U)&(I - n*n))
	  + lambdaf*tr(shearGradU&(I - n*n))*n
	  + muf*(shearGradU&n)
	  )
	);
   }
 else if(divDSigmaExpMethod == "expLaplacian")
   {
     divDSigmaExp =
       - fvc::laplacian(mu + lambda, U, "laplacian(DU,U)")
       + fvc::div
       (
	mu*gradDU.T()
	+ lambda*(I*tr(gradDU)),
	"div(sigma)"
	);
   }
 else
   {
     FatalError << "divDSigmaExp method " << divDSigmaExpMethod << " not found!" << endl;
   }
