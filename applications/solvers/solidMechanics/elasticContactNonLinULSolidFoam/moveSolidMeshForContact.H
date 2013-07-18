{
  // Create point interpolation
  volPointInterpolation pointInterpolation(mesh);

  // Calculate mesh points displacement
  pointVectorField pointDU = pointInterpolation.interpolate(DU);

  vectorField pointDUI = pointDU.internalField();

  // Move mesh
  vectorField newPoints = mesh.allPoints();

  forAll (pointDUI, pointI)
    {
      newPoints[pointI] += pointDUI[pointI];
    }

# include "correctGlobalFaceZoneMesh.H"

  twoDPointCorrector twoDCorrector(mesh);
  twoDCorrector.correctPoints(newPoints);

  mesh.movePoints(newPoints);
  //  pMesh.movePoints(newPoints);
  mesh.V00();
  mesh.moving(false);
}
