
// Create solid solver
autoPtr<solidSolver> stress = solidSolver::New(mesh);

// label nWedgePatches = 0;
// vector wedgeDirVec = vector::zero;
// forAll(mesh.boundaryMesh(), patchI)
// {
//     if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
//     {
//         const wedgePolyPatch& wpp = refCast<const wedgePolyPatch>
//             (
//                 mesh.boundaryMesh()[patchI]
//             );

//         nWedgePatches++;
//         wedgeDirVec += cmptMag(wpp.centreNormal());
//     }
// }

// reduce(nWedgePatches, maxOp<label>());

// if (nWedgePatches)
// {
//     Info<< nl << "Axisymmetric case: disabling the 3rd direction"
//         << nl << endl;

//     // We will const_cast as it is currently the tidiest way, until polyMesh is
//     // modified or gives write access to solutionD
//     Vector<label>& solD = const_cast<Vector<label>&>(mesh.solutionD());

//     reduce(wedgeDirVec, sumOp<vector>());

//     wedgeDirVec /= mag(wedgeDirVec);

//     for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
//     {
//         if (wedgeDirVec[cmpt] > 1e-6)
//         {
//             solD[cmpt] = -1;
//         }
//         else
//         {
//             solD[cmpt] = 1;
//         }
//     }
// }
