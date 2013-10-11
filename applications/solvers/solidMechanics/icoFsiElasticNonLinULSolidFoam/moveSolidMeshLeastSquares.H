//--------------------------------------------------//
//- move mesh
//--------------------------------------------------//
if(min(J.internalField()) > 0)
  {
    Info << "Moving mesh using least squares interpolation" << endl;

    leastSquaresVolPointInterpolation pointInterpolation(stressMesh);
  
    // Create point mesh
    pointMesh pMesh(stressMesh);
    
    wordList types
      (
       pMesh.boundary().size(),
       calculatedFvPatchVectorField::typeName
       );
    
    pointVectorField pointDU
      (
       IOobject
       (
        "pointDU",
        runTime.timeName(),
        stressMesh
        ),
       pMesh,
       dimensionedVector("zero", dimLength, vector::zero),
       types
       );
    
    pointInterpolation.interpolate(DU, pointDU);

    const vectorField& pointDUI =
      pointDU.internalField();

    //- Move mesh
    vectorField newPoints = stressMesh.allPoints();

    forAll (pointDUI, pointI)
      {
        newPoints[pointI] += pointDUI[pointI];
      }

    twoDPointCorrector twoDCorrector(stressMesh);
    twoDCorrector.correctPoints(newPoints);
    stressMesh.movePoints(newPoints);
    stressMesh.V00();
    stressMesh.moving(false);
  }
  else
    {
      FatalErrorIn(args.executable())
	<< "Negative Jacobian"
	<< exit(FatalError);
    }
