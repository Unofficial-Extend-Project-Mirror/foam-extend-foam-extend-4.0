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

    volScalarField pressure
      (
       IOobject
       (
        "pressure",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
       tr(sigma)/3.0
       );
    //- boundary surface pressure
    forAll(pressure.boundaryField(), patchi)
      {
        const vectorField& nb = n.boundaryField()[patchi];
        pressure.boundaryField()[patchi] =
          -(nb & ( nb & sigma.boundaryField()[patchi] ));
      }

    //- contact slave penetration
    #     include "moveSolidMesh.H"
    pointMesh pMesh(mesh);
    pointScalarField cGapPoints
      (
       IOobject
       (
        "pointContactGap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
       pMesh,
       dimensionedScalar("scalar", dimLength, 0.0),
       "calculated"
       );
    contact.contactGapPoints(cGapPoints);
    cGapPoints.write();
    mesh.movePoints(oldMeshPoints);  

    runTime.write();
  }
