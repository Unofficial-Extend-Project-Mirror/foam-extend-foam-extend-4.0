//--------------------------------------------------//
//- move mesh
//--------------------------------------------------//
// if(min(J.internalField()) > 0)
{
    Info << "Moving mesh using least squares interpolation" << endl;

    leastSquaresVolPointInterpolation pointInterpolation(mesh);

    // Create point mesh
    pointMesh pMesh(mesh);

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

    const vectorField& pointDUI = pointDU.internalField();

    //- Move mesh
    vectorField newPoints = mesh.allPoints();

    forAll (pointDUI, pointI)
    {
        newPoints[pointI] += pointDUI[pointI];
    }

    twoDPointCorrector twoDCorrector(mesh);
    twoDCorrector.correctPoints(newPoints);
    mesh.movePoints(newPoints);
    mesh.V00();
    mesh.moving(false);
}
  // else
  //   {
  //     FatalErrorIn(args.executable())
  //     << "Negative Jacobian"
  //     << exit(FatalError);
  //   }
