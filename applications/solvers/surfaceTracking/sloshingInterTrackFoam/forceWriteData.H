    OFstream file
    (
        runTime.timePath()/"curvature.dat"
    );

    const vectorField& C = interface.aMesh().areaCentres().internalField();
    const scalarField& K = interface.aMesh().faceCurvatures().internalField();

//     const scalarField divUs = fac::div(Us)().internalField();

// const vectorField& nP = interface.aMesh().pointAreaNormals();
// const vectorField& nF = interface.aMesh().faceAreaNormals();

//     primitivePatchInterpolationZT 
//     patchInterpolator(mesh.boundaryMesh()[movingPatchID]);

//     vectorField nFAvr =
//         patchInterpolator.pointToFaceInterpolate(nP);

//     nFAvr /= mag(nFAvr);

    forAll(C, faceI)
    {
        file<< C[faceI].x() << ' ' 
            << C[faceI].y() << ' '
            << K[faceI] << endl;
//             << phi.boundaryField()[movingPatchID][faceI] << ' '
//             << (180.0/M_PI)*std::acos((nF[faceI])&(nFAvr[faceI])) << endl;
    }
