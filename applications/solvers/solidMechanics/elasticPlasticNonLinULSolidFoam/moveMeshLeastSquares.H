//--------------------------------------------------//
//- move mesh
//--------------------------------------------------//
{
    //Info << "Moving mesh using least squares interpolation" << endl;

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

    const vectorField& pointDUI = pointDU.internalField();

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

    // note: allPoints will have more points than pointDU
    // if there are globalFaceZones
    forAll (pointDUI, pointI)
    {
        newPoints[pointI] += pointDUI[pointI];
    }

    // Correct symmetryPlane points

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<symmetryPolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            vector avgN =
                gAverage(mesh.boundaryMesh()[patchI].pointNormals());

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

//     // now we make sure processor patches are exactly the same
//     newPoints.correctBoundaryConditions();

    twoDPointCorrector twoDCorrector(mesh);
    twoDCorrector.correctPoints(newPoints);
    mesh.movePoints(newPoints);
    mesh.V00();
    mesh.moving(false);

    // Update n
    n = mesh.Sf()/mesh.magSf();
}
