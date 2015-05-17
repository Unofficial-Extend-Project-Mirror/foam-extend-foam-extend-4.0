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

    // Correct symmetryPlane points

    forAll(stressMesh.boundaryMesh(), patchI)
    {
        if (isA<symmetryPolyPatch>(stressMesh.boundaryMesh()[patchI]))
        {
            const labelList& meshPoints =
                stressMesh.boundaryMesh()[patchI].meshPoints();

            vector avgN =
                gAverage(stressMesh.boundaryMesh()[patchI].pointNormals());

            vector i(1, 0, 0);
            vector j(0, 1, 0);
            vector k(0, 0, 1);

            if (mag(avgN&i) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    newPoints[meshPoints[pI]].x() = 0;
                }
            }
            else if (mag(avgN&j) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    newPoints[meshPoints[pI]].y() = 0;
                }
            }
            else if (mag(avgN&k) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    newPoints[meshPoints[pI]].z() = 0;
                }
            }
        }
    }

#   include "calcUnusedNewPoints.H"

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
