// Calculate third order dissipative term as outlined by Demirdzic
// This removes unphysical high frequency oscillations that may
// appear in the solution. This term goes to zero on convergence
// assuming a linear/quadratic displacement field, or the error
// is second order

// volVectorField divThirdOrderTerm
// (
//     IOobject
//     (
//         "divThirdOrderTerm",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedVector("zero", dimForce/dimVolume, vector::zero)
// );

// average gradU of neighbouring cell centres
// interpolation scheme should be midPoint
surfaceTensorField averageGradU("averageGradU", fvc::interpolate(gradU, "averageGradU"));

// average face gradU extrapolated from neighbouring cell centres
surfaceTensorField extrapGradU
(
    IOobject
    (
        "extrapGradU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("zero", dimless, tensor::zero)
);

volVectorField gradGradUcompXX = fvc::grad(gradU.component(tensor::XX), "gradGradU");
volVectorField gradGradUcompXY = fvc::grad(gradU.component(tensor::XY), "gradGradU");
volVectorField gradGradUcompXZ = fvc::grad(gradU.component(tensor::XZ), "gradGradU");
volVectorField gradGradUcompYX = fvc::grad(gradU.component(tensor::YX), "gradGradU");
volVectorField gradGradUcompYY = fvc::grad(gradU.component(tensor::YY), "gradGradU");
volVectorField gradGradUcompYZ = fvc::grad(gradU.component(tensor::YZ), "gradGradU");
volVectorField gradGradUcompZX = fvc::grad(gradU.component(tensor::ZX), "gradGradU");
volVectorField gradGradUcompZY = fvc::grad(gradU.component(tensor::ZY), "gradGradU");
volVectorField gradGradUcompZZ = fvc::grad(gradU.component(tensor::ZZ), "gradGradU");

// third order tensor components
volScalarField gradGradUXXX = gradGradUcompXX.component(vector::X);
volScalarField gradGradUXXY = gradGradUcompXY.component(vector::X);
volScalarField gradGradUXXZ = gradGradUcompXZ.component(vector::X);

volScalarField gradGradUXYX = gradGradUcompYX.component(vector::X);
volScalarField gradGradUXYY = gradGradUcompYY.component(vector::X);
volScalarField gradGradUXYZ = gradGradUcompYZ.component(vector::X);

volScalarField gradGradUXZX = gradGradUcompZX.component(vector::X);
volScalarField gradGradUXZY = gradGradUcompZY.component(vector::X);
volScalarField gradGradUXZZ = gradGradUcompZZ.component(vector::X);

volScalarField gradGradUYXX = gradGradUcompXX.component(vector::Y);
volScalarField gradGradUYXY = gradGradUcompXY.component(vector::Y);
volScalarField gradGradUYXZ = gradGradUcompXZ.component(vector::Y);

volScalarField gradGradUYYX = gradGradUcompYX.component(vector::Y);
volScalarField gradGradUYYY = gradGradUcompYY.component(vector::Y);
volScalarField gradGradUYYZ = gradGradUcompYZ.component(vector::Y);

volScalarField gradGradUYZX = gradGradUcompZX.component(vector::Y);
volScalarField gradGradUYZY = gradGradUcompZY.component(vector::Y);
volScalarField gradGradUYZZ = gradGradUcompZZ.component(vector::Y);

volScalarField gradGradUZXX = gradGradUcompXX.component(vector::Z);
volScalarField gradGradUZXY = gradGradUcompXY.component(vector::Z);
volScalarField gradGradUZXZ = gradGradUcompXZ.component(vector::Z);

volScalarField gradGradUZYX = gradGradUcompYX.component(vector::Z);
volScalarField gradGradUZYY = gradGradUcompYY.component(vector::Z);
volScalarField gradGradUZYZ = gradGradUcompYZ.component(vector::Z);

volScalarField gradGradUZZX = gradGradUcompZX.component(vector::Z);
volScalarField gradGradUZZY = gradGradUcompZY.component(vector::Z);
volScalarField gradGradUZZZ = gradGradUcompZZ.component(vector::Z);

forAll(extrapGradU.internalField(), facei)
{
    const label own = mesh.owner()[facei];
    const label nei = mesh.neighbour()[facei];
    const vector deltaOwn = mesh.Cf()[facei] - mesh.C()[own];
    const vector deltaNei = mesh.Cf()[facei] - mesh.C()[nei];
    const tensor& gradUOwn = gradU.internalField()[own];
    const tensor& gradUNei = gradU.internalField()[nei];

    // as there is there is no thirdOrderTensor class, we will
    // calculate (deltaOwn&gradGradUOwn) out manually
    // tensor deltaOwnDotgradGradUOwn = tensor::zero;
    // tensor deltaNeiDotgradGradUNei = tensor::zero;

    // deltaOwnDotgradGradUOwn[tensor::XX] = deltaOwn & gradGradUcompXX.internalField()[own];
    // deltaNeiDotgradGradUNei[tensor::XX] = deltaNei & gradGradUcompXX.internalField()[nei];
    // deltaOwnDotgradGradUOwn[tensor::XY] = deltaOwn & gradGradUcompXY.internalField()[own];
    // deltaNeiDotgradGradUNei[tensor::XY] = deltaNei & gradGradUcompXY.internalField()[nei];
    // deltaOwnDotgradGradUOwn[tensor::XZ] = deltaOwn & gradGradUcompXZ.internalField()[own];
    // deltaNeiDotgradGradUNei[tensor::XZ] = deltaNei & gradGradUcompXZ.internalField()[nei];

    // deltaOwnDotgradGradUOwn[tensor::YX] = deltaOwn & gradGradUcompYX.internalField()[own];
    // deltaNeiDotgradGradUNei[tensor::YX] = deltaNei & gradGradUcompYX.internalField()[nei];
    // deltaOwnDotgradGradUOwn[tensor::YY] = deltaOwn & gradGradUcompYY.internalField()[own];
    // deltaNeiDotgradGradUNei[tensor::YY] = deltaNei & gradGradUcompYY.internalField()[nei];
    // deltaOwnDotgradGradUOwn[tensor::YZ] = deltaOwn & gradGradUcompYZ.internalField()[own];
    // deltaNeiDotgradGradUNei[tensor::YZ] = deltaNei & gradGradUcompYZ.internalField()[nei];

    // deltaOwnDotgradGradUOwn[tensor::ZX] = deltaOwn & gradGradUcompZX.internalField()[own];
    // deltaNeiDotgradGradUNei[tensor::ZX] = deltaNei & gradGradUcompZX.internalField()[nei];
    // deltaOwnDotgradGradUOwn[tensor::ZY] = deltaOwn & gradGradUcompZY.internalField()[own];
    // deltaNeiDotgradGradUNei[tensor::ZY] = deltaNei & gradGradUcompZY.internalField()[nei];
    // deltaOwnDotgradGradUOwn[tensor::ZZ] = deltaOwn & gradGradUcompZZ.internalField()[own];
    // deltaNeiDotgradGradUNei[tensor::ZZ] = deltaNei & gradGradUcompZZ.internalField()[nei];

    scalar gradGradUXXXOwn = gradGradUXXX.internalField()[own];
    scalar gradGradUXXYOwn = gradGradUXXY.internalField()[own];
    scalar gradGradUXXZOwn = gradGradUXXZ.internalField()[own];

    scalar gradGradUXYXOwn = gradGradUXYX.internalField()[own];
    scalar gradGradUXYYOwn = gradGradUXYY.internalField()[own];
    scalar gradGradUXYZOwn = gradGradUXYZ.internalField()[own];

    scalar gradGradUXZXOwn = gradGradUXZX.internalField()[own];
    scalar gradGradUXZYOwn = gradGradUXZY.internalField()[own];
    scalar gradGradUXZZOwn = gradGradUXZZ.internalField()[own];


    scalar gradGradUYXXOwn = gradGradUYXX.internalField()[own];
    scalar gradGradUYXYOwn = gradGradUYXY.internalField()[own];
    scalar gradGradUYXZOwn = gradGradUYXZ.internalField()[own];

    scalar gradGradUYYXOwn = gradGradUYYX.internalField()[own];
    scalar gradGradUYYYOwn = gradGradUYYY.internalField()[own];
    scalar gradGradUYYZOwn = gradGradUYYZ.internalField()[own];

    scalar gradGradUYZXOwn = gradGradUYZX.internalField()[own];
    scalar gradGradUYZYOwn = gradGradUYZY.internalField()[own];
    scalar gradGradUYZZOwn = gradGradUYZZ.internalField()[own];


    scalar gradGradUZXXOwn = gradGradUZXX.internalField()[own];
    scalar gradGradUZXYOwn = gradGradUZXY.internalField()[own];
    scalar gradGradUZXZOwn = gradGradUZXZ.internalField()[own];

  scalar gradGradUZYXOwn = gradGradUZYX.internalField()[own];
  scalar gradGradUZYYOwn = gradGradUZYY.internalField()[own];
  scalar gradGradUZYZOwn = gradGradUZYZ.internalField()[own];

  scalar gradGradUZZXOwn = gradGradUZZX.internalField()[own];
  scalar gradGradUZZYOwn = gradGradUZZY.internalField()[own];
  scalar gradGradUZZZOwn = gradGradUZZZ.internalField()[own];


    // nei
    scalar gradGradUXXXNei = gradGradUXXX.internalField()[nei];
    scalar gradGradUXXYNei = gradGradUXXY.internalField()[nei];
    scalar gradGradUXXZNei = gradGradUXXZ.internalField()[nei];

    scalar gradGradUXYXNei = gradGradUXYX.internalField()[nei];
    scalar gradGradUXYYNei = gradGradUXYY.internalField()[nei];
    scalar gradGradUXYZNei = gradGradUXYZ.internalField()[nei];

    scalar gradGradUXZXNei = gradGradUXZX.internalField()[nei];
    scalar gradGradUXZYNei = gradGradUXZY.internalField()[nei];
    scalar gradGradUXZZNei = gradGradUXZZ.internalField()[nei];


    scalar gradGradUYXXNei = gradGradUYXX.internalField()[nei];
    scalar gradGradUYXYNei = gradGradUYXY.internalField()[nei];
    scalar gradGradUYXZNei = gradGradUYXZ.internalField()[nei];

    scalar gradGradUYYXNei = gradGradUYYX.internalField()[nei];
    scalar gradGradUYYYNei = gradGradUYYY.internalField()[nei];
    scalar gradGradUYYZNei = gradGradUYYZ.internalField()[nei];

    scalar gradGradUYZXNei = gradGradUYZX.internalField()[nei];
    scalar gradGradUYZYNei = gradGradUYZY.internalField()[nei];
    scalar gradGradUYZZNei = gradGradUYZZ.internalField()[nei];


    scalar gradGradUZXXNei = gradGradUZXX.internalField()[nei];
    scalar gradGradUZXYNei = gradGradUZXY.internalField()[nei];
    scalar gradGradUZXZNei = gradGradUZXZ.internalField()[nei];

    scalar gradGradUZYXNei = gradGradUZYX.internalField()[nei];
    scalar gradGradUZYYNei = gradGradUZYY.internalField()[nei];
    scalar gradGradUZYZNei = gradGradUZYZ.internalField()[nei];

    scalar gradGradUZZXNei = gradGradUZZX.internalField()[nei];
    scalar gradGradUZZYNei = gradGradUZZY.internalField()[nei];
    scalar gradGradUZZZNei = gradGradUZZZ.internalField()[nei];

    tensor deltaOwnDotgradGradUOwn =tensor
    (
        deltaOwn.x()*gradGradUXXXOwn + deltaOwn.y()*gradGradUYXXOwn + deltaOwn.z()*gradGradUZXXOwn,
        deltaOwn.x()*gradGradUXXYOwn + deltaOwn.y()*gradGradUYXYOwn + deltaOwn.z()*gradGradUZXYOwn,
        deltaOwn.x()*gradGradUXXZOwn + deltaOwn.y()*gradGradUYXZOwn + deltaOwn.z()*gradGradUZXZOwn,

        deltaOwn.x()*gradGradUXYXOwn + deltaOwn.y()*gradGradUYYXOwn + deltaOwn.z()*gradGradUZYXOwn,
        deltaOwn.x()*gradGradUXYYOwn + deltaOwn.y()*gradGradUYYYOwn + deltaOwn.z()*gradGradUZYYOwn,
        deltaOwn.x()*gradGradUXYZOwn + deltaOwn.y()*gradGradUYYZOwn + deltaOwn.z()*gradGradUZYZOwn,

        deltaOwn.x()*gradGradUXZXOwn + deltaOwn.y()*gradGradUYZXOwn + deltaOwn.z()*gradGradUZZXOwn,
        deltaOwn.x()*gradGradUXZYOwn + deltaOwn.y()*gradGradUYZYOwn + deltaOwn.z()*gradGradUZZYOwn,
        deltaOwn.x()*gradGradUXZZOwn + deltaOwn.y()*gradGradUYZZOwn + deltaOwn.z()*gradGradUZZZOwn
    );

    tensor deltaNeiDotgradGradUNei = tensor
    (
        deltaNei.x()*gradGradUXXXNei + deltaNei.y()*gradGradUYXXNei + deltaNei.z()*gradGradUZXXNei,
        deltaNei.x()*gradGradUXXYNei + deltaNei.y()*gradGradUYXYNei + deltaNei.z()*gradGradUZXYNei,
        deltaNei.x()*gradGradUXXZNei + deltaNei.y()*gradGradUYXZNei + deltaNei.z()*gradGradUZXZNei,

        deltaNei.x()*gradGradUXYXNei + deltaNei.y()*gradGradUYYXNei + deltaNei.z()*gradGradUZYXNei,
        deltaNei.x()*gradGradUXYYNei + deltaNei.y()*gradGradUYYYNei + deltaNei.z()*gradGradUZYYNei,
        deltaNei.x()*gradGradUXYZNei + deltaNei.y()*gradGradUYYZNei + deltaNei.z()*gradGradUZYZNei,

        deltaNei.x()*gradGradUXZXNei + deltaNei.y()*gradGradUYZXNei + deltaNei.z()*gradGradUZZXNei,
        deltaNei.x()*gradGradUXZYNei + deltaNei.y()*gradGradUYZYNei + deltaNei.z()*gradGradUZZYNei,
        deltaNei.x()*gradGradUXZZNei + deltaNei.y()*gradGradUYZZNei + deltaNei.z()*gradGradUZZZNei
    );


    // get average of extrapolated values
    tensor extrapFaceGrad =
        0.5*
        (
            gradUOwn + (deltaOwnDotgradGradUOwn)
            +
            gradUNei + (deltaNeiDotgradGradUNei)
        );

    extrapGradU.internalField()[facei] = extrapFaceGrad;

}

// correction is zero on boundary
forAll(extrapGradU.boundaryField(), patchi)
{
    extrapGradU.boundaryField()[patchi] =
        averageGradU.boundaryField()[patchi];
}

// calculate thirdOrderTerm
volVectorField divThirdOrderTerm
(
    "thirdOrderTerm",
    fvc::div
    (
        (2*muf+lambdaf)*mesh.Sf() & (extrapGradU - averageGradU)
    )
);

// if(runTime.outputTime())
// {
//     divThirdOrderTerm.write();
//     averageGradU.write();
//     extrapGradU.write();
// }
