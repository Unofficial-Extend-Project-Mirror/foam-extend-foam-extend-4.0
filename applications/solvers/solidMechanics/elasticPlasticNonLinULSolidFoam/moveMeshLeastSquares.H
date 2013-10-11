//--------------------------------------------------//
//- move mesh
//--------------------------------------------------//
  {
    Info << "Moving mesh using least squares interpolation" << endl;
    
    leastSquaresVolPointInterpolation pointInterpolation(mesh);
  
    // Create point mesh
    pointMesh pMesh(mesh);
    //pointMesh pMesh(mesh);
    
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
        mesh
        ),
       pMesh,
       dimensionedVector("zero", dimLength, vector::zero),
       types
       );

    pointInterpolation.interpolate(DU, pointDU);

    //pointDU.write();

    const vectorField& pointDUI =
      pointDU.internalField();

    //- Move mesh
    //vectorField newPoints = mesh.allPoints();
    pointVectorField newPoints
      (
       IOobject
       (
        "newPoints",
        runTime.timeName(),
        mesh
        ),
       pMesh,
       dimensionedVector("zero", dimLength, vector::zero)
       //mesh.allPoints()
       );
    newPoints.internalField() = mesh.allPoints();
    //const vectorField oldPoints = mesh.allPoints();

    // note: allPoints will have more points than pointDU
    // if there are globalFaceZones
    forAll (pointDUI, pointI)
      {
        newPoints[pointI] += pointDUI[pointI];
      }

    // now we make sure processor patches are exactly the same
    newPoints.correctBoundaryConditions();
    
    twoDPointCorrector twoDCorrector(mesh);
    twoDCorrector.correctPoints(newPoints);
    mesh.movePoints(newPoints);
    mesh.V00();
    mesh.moving(false);

  }
