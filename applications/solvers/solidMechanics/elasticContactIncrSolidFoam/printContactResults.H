if (runTime.outputTime())
  {
    // FAILS IN PARALLEL - FIX
    //    Info << "Print contact area" << endl;
    //volScalarField ca = contact.contactArea();
    //ca.write();
    
    //-------------------------------------------------------------//
    // I couldn't get tmp to return the pointScalarField correctly //
    // so I had to make the pointScalarField here and pass it to   //
    // contactGapPoints and pointContactForce to populate          //
    //-------------------------------------------------------------//
    //This is the point distance for each contact vertex
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
	  
    
    //- This is the point distance for each contact vertex
    pointVectorField cPointForce 
      (
       IOobject
       (
	"pointContactForce",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
       pMesh,
       dimensionedVector("vector", dimForce, vector::zero),
       "calculated"
       );
    contact.contactPointForce(cPointForce);
    cPointForce.write();
    
    //- this is the actual (sigma&n)&n) on the contact patches
    //- SHOULD THIS BE A REF TO A TMP...?
    volScalarField cPressure = contact.contactPressure();
    cPressure.write();
  }
