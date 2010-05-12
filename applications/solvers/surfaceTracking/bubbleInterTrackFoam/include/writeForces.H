if(Pstream::master())
{
//     const scalarField& magSf =
//         mesh.magSf().boundaryField()[interface.aPatchID()];

//     vectorField n =
//         mesh.Sf().boundaryField()[interface.aPatchID()]/magSf;

//     const scalarField& P = p.boundaryField()[interface.aPatchID()];
 
//     vectorField pressureForces = magSf*P*n;

//     vectorField snGradU = 
//         U.boundaryField()[interface.aPatchID()].snGrad();

//     vectorField viscousForces = 
//       - interface.muFluidA().value()*magSf*
//         (
//             snGradU 
//           + ((n*n)&snGradU) 
//           + (fac::grad(interface.Us())&interface.aMesh().faceAreaNormals())
//             ().internalField()
//         );

    vector totalPressureForce = -interface.totalPressureForce();

    vector totalViscousForce = -interface.totalViscousForce();

    vector totalForce = totalPressureForce + totalViscousForce;


    vector dragDir;
    scalar Uref;

    if(mag(UF.value()) > SMALL)
    {
        dragDir = - UF.value()/mag(UF.value());
        Uref = mag(UF.value());
    }
    else
    {
        dragDir = interface.g().value()/mag(interface.g().value());
        Uref = SMALL;
    }


    scalar Aref = 
        sum
        (
            (dragDir&mesh.Sf().boundaryField()[interface.aPatchID()])*
            pos(dragDir&mesh.Sf().boundaryField()[interface.aPatchID()])
        );


    scalar drag =
        (dragDir&totalForce)/
        (0.5*interface.rhoFluidA().value()*sqr(Uref)*Aref);


    scalar lift =
        mag((I - dragDir*dragDir)&totalForce)/
        (0.5*interface.rhoFluidA().value()*sqr(Uref)*Aref);


    forcesFile << runTime.value() << tab;

    forcesFile << drag << tab << lift << endl;
}
