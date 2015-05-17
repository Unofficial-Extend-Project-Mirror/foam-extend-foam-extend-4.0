if (runTime.outputTime())
{
    //C.write();

    volScalarField epsilonEq
    (
        IOobject
        (
            "epsilonEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((2.0/3.0)*magSqr(dev(epsilon)))
    );

    Info<< "Max epsilonEq = " << max(epsilonEq).value() << endl;

    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma)))
    );

    Info<< "Max sigmaEq = " << max(sigmaEq).value() << endl;

// volVectorField traction
// (
//     IOobject
//     (
//         "traction",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedVector("zero", dimForce/dimArea, vector::zero)
// );
// forAll(mesh.boundary(), patchi)
// {
//     traction.boundaryField()[patchi] =
//         n.boundaryField()[patchi] & sigma.boundaryField()[patchi];
// }

//     //- patch forces
//     forAll(mesh.boundary(), patchi)
//     {
//         Info << "Patch " << mesh.boundary()[patchi].name() << endl;
//         vectorField totalForce = mesh.Sf().boundaryField()[patchi] & sigma.boundaryField()[patchi];
//         vector force = sum( totalForce );
//         Info << "\ttotal force is " << force << " N" << endl;
//     tensorField F = I + gradDU.boundaryField()[patchi];
//         tensorField Finv = inv(F);
//         scalar normalForce = sum( n.boundaryField()[patchi] & totalForce );
//         Info << "\tnormal force is " << normalForce << " N" << endl;
//         scalar shearForce = mag(sum( (I - sqr(n.boundaryField()[patchi])) & totalForce ));
//         Info << "\tshear force is " << shearForce << " N" << endl;
//     }

    runTime.write();
}
