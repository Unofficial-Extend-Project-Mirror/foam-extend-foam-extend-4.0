Info << " Reading stress mesh" << endl;

// Info << runTime.timePath() << endl;

fvMesh solidMesh
(
    IOobject
    (
        "solid",
//         runTime.constant(),
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);


// // Correct symmetryPlane points

// vectorField& newPoints = const_cast<vectorField&>(solidMesh.points());

// forAll(solidMesh.boundaryMesh(), patchI)
// {
//     if (isA<symmetryPolyPatch>(solidMesh.boundaryMesh()[patchI]))
//     {
//         const labelList& meshPoints =
//         solidMesh.boundaryMesh()[patchI].meshPoints();

//     vector avgN =
//         gAverage(solidMesh.boundaryMesh()[patchI].pointNormals());

//     vector i(1, 0, 0);
//     vector j(0, 1, 0);
//     vector k(0, 0, 1);

//     if (mag(avgN&i) > 0.95)
//     {
//         forAll(meshPoints, pI)
//         {
//             newPoints[meshPoints[pI]].x() = 0;
//         }
//     }
//     else if (mag(avgN&j) > 0.95)
//     {
//         forAll(meshPoints, pI)
//         {
//             newPoints[meshPoints[pI]].y() = 0;
//         }
//     }
//     else if (mag(avgN&k) > 0.95)
//     {
//         forAll(meshPoints, pI)
//         {
//             newPoints[meshPoints[pI]].z() = 0;
//         }
//     }

//     avgN = gAverage(solidMesh.boundaryMesh()[patchI].pointNormals());

//     Pout << "avgN = " << avgN << endl;
//     }
// }

// solidMesh.movePoints(newPoints);
// solidMesh.V00();
// solidMesh.moving(false);
