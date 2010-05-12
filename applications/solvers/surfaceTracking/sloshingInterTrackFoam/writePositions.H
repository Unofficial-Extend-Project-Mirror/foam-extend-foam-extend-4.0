
    positionsFile << runTime.value() << tab;

    const pointField& P = 
        mesh.boundaryMesh()[interface.aPatchID()].localPoints();

//     scalar ePot = 
//         mag
//         (
//             rhoLiquid.value()*
//             sum((g.value()&mesh.C().internalField())*mesh.V())
//         );

//     scalar eKin = 
//         0.5*rhoLiquid.value()*sum(magSqr(U)().internalField()*mesh.V());


//     tensorField gradU = fvc::grad(U, "test")().internalField();

//     dissipation += 
//         2.0*muLiquid.value()
//         *sum((symm(gradU) && gradU)*mesh.V())
//         *runTime.deltaT().value();


//     fvScalarMatrix eTotEqn
//     (
//         fvm::ddt(eTot)
//       + fvm::div(phi, eTot)
//       + dissipation
//     );

//     eTotEqn.solve();

    positionsFile.precision(12);

positionsFile << P[leftPointLabel].y() << tab 
    << P[rightPointLabel].y() << endl;
//         << ePot << tab
//         << eKin << tab
//         << ePot + eKin << tab
//         << ePot + eKin + dissipation << endl;

