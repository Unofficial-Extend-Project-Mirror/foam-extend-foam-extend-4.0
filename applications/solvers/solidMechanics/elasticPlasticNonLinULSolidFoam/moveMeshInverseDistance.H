//--------------------------------------------------//
//- move mesh
//--------------------------------------------------//
// if(min(J.internalField()) > 0)
{
    Info << "Move solid mesh using inverse distance interpolation" << endl;

    // Create point mesh
    pointMesh pMesh(mesh);

    // Create point interpolation
    volPointInterpolation pointInterpolation(mesh);

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

    // Calculate mesh points displacement
    pointInterpolation.interpolate(DU, pointDU);

    //- correct edge interpolation
    //- this is the stuff from edgeCorrectedVolPointInterpolation but
    //- that class no longer works
    # include "performEdgeCorrectedVolPointInterpolation.H"

    const vectorField& pointDUI =
        pointDU.internalField();

    //- see the effect of correctBCs

    // Move mesh
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
//   else
//     {
//       FatalErrorIn(args.executable())
//     << "Negative Jacobian"
//     << exit(FatalError);
//     }
