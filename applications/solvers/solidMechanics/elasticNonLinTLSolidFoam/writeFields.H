if (runTime.outputTime())
  {
    volScalarField epsilonEq
      (
       IOobject
       (
	"epsilonEq",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
       sqrt((2.0/3.0)*magSqr(dev(epsilon)))
       );

    Info<< "Max epsilonEq = " << max(epsilonEq).value()
	<< endl;

    volScalarField sigmaEq
      (
       IOobject
       (
	"sigmaEq",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
       sqrt((3.0/2.0)*magSqr(dev(sigma)))
       );

    Info<< "Max sigmaEq = " << max(sigmaEq).value()
	<< endl;

    //- Calculate Cauchy stress
    volTensorField F = I + gradU;
    volScalarField J = det(F);

    //- update density
    rho = rho/J;

    volSymmTensorField sigmaCauchy
      (
       IOobject
       (
	"sigmaCauchy",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
       (1/J) * symm(F.T() & sigma & F)
       );

    //- Cauchy von Mises stress
    volScalarField sigmaCauchyEq
      (
       IOobject
       (
	"sigmaCauchyEq",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
       sqrt((3.0/2.0)*magSqr(dev(sigmaCauchy)))
       );

    Info<< "Max sigmaCauchyEq = " << max(sigmaCauchyEq).value()
	<< endl;

    //- write boundary forces
    //- integrate (sigma2PK & F) over reference area
    //- which is equivalent to integrating sigmaCauchy
    //- over the deformed area
    Info << nl;
    forAll(mesh.boundary(), patchi)
      {
	Info << "Patch " << mesh.boundary()[patchi].name() << endl;
	tensorField F = I + gradU.boundaryField()[patchi];
	vectorField totalForce = mesh.Sf().boundaryField()[patchi] & (sigma.boundaryField()[patchi] & F);

	vector force = sum( totalForce );
	Info << "\ttotal force is " << force << " N" << endl;

	tensorField Finv = inv(F);
	vectorField nCurrent = Finv & n.boundaryField()[patchi];
	nCurrent /= mag(nCurrent);
	scalar normalForce = sum( nCurrent & totalForce );
	Info << "\tnormal force is " << normalForce << " N" << endl;
	scalar shearForce = mag(sum( (I - sqr(nCurrent)) & totalForce ));
	Info << "\tshear force is " << shearForce << " N" << endl << endl;;
      }


    //- move mesh for visualisation and move it back after writing
    vectorField oldPoints = mesh.allPoints();
    #include "moveMeshLeastSquares.H"

    runTime.write();

    //- move mesh back
    mesh.movePoints(oldPoints);
  }
