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

// average gradDU of neighbouring cell centres
// interpolation scheme should be midPoint
surfaceTensorField averageGradDU("averageGradDU", fvc::interpolate(gradDU, "averageGradDU"));

// average face gradDU extrapolated from neighbouring cell centres
surfaceTensorField extrapGradDU
(
    IOobject
    (
        "extrapGradDU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("zero", dimless, tensor::zero)
);

volVectorField gradGradDUcompXX = fvc::grad(gradDU.component(tensor::XX), "gradGradDU");
volVectorField gradGradDUcompXY = fvc::grad(gradDU.component(tensor::XY), "gradGradDU");
volVectorField gradGradDUcompXZ = fvc::grad(gradDU.component(tensor::XZ), "gradGradDU");
volVectorField gradGradDUcompYX = fvc::grad(gradDU.component(tensor::YX), "gradGradDU");
volVectorField gradGradDUcompYY = fvc::grad(gradDU.component(tensor::YY), "gradGradDU");
volVectorField gradGradDUcompYZ = fvc::grad(gradDU.component(tensor::YZ), "gradGradDU");
volVectorField gradGradDUcompZX = fvc::grad(gradDU.component(tensor::ZX), "gradGradDU");
volVectorField gradGradDUcompZY = fvc::grad(gradDU.component(tensor::ZY), "gradGradDU");
volVectorField gradGradDUcompZZ = fvc::grad(gradDU.component(tensor::ZZ), "gradGradDU");

// third order tensor components
volScalarField gradGradDUXXX = gradGradDUcompXX.component(vector::X);
volScalarField gradGradDUXXY = gradGradDUcompXY.component(vector::X);
volScalarField gradGradDUXXZ = gradGradDUcompXZ.component(vector::X);

volScalarField gradGradDUXYX = gradGradDUcompYX.component(vector::X);
volScalarField gradGradDUXYY = gradGradDUcompYY.component(vector::X);
volScalarField gradGradDUXYZ = gradGradDUcompYZ.component(vector::X);

volScalarField gradGradDUXZX = gradGradDUcompZX.component(vector::X);
volScalarField gradGradDUXZY = gradGradDUcompZY.component(vector::X);
volScalarField gradGradDUXZZ = gradGradDUcompZZ.component(vector::X);

volScalarField gradGradDUYXX = gradGradDUcompXX.component(vector::Y);
volScalarField gradGradDUYXY = gradGradDUcompXY.component(vector::Y);
volScalarField gradGradDUYXZ = gradGradDUcompXZ.component(vector::Y);

volScalarField gradGradDUYYX = gradGradDUcompYX.component(vector::Y);
volScalarField gradGradDUYYY = gradGradDUcompYY.component(vector::Y);
volScalarField gradGradDUYYZ = gradGradDUcompYZ.component(vector::Y);

volScalarField gradGradDUYZX = gradGradDUcompZX.component(vector::Y);
volScalarField gradGradDUYZY = gradGradDUcompZY.component(vector::Y);
volScalarField gradGradDUYZZ = gradGradDUcompZZ.component(vector::Y);

volScalarField gradGradDUZXX = gradGradDUcompXX.component(vector::Z);
volScalarField gradGradDUZXY = gradGradDUcompXY.component(vector::Z);
volScalarField gradGradDUZXZ = gradGradDUcompXZ.component(vector::Z);

volScalarField gradGradDUZYX = gradGradDUcompYX.component(vector::Z);
volScalarField gradGradDUZYY = gradGradDUcompYY.component(vector::Z);
volScalarField gradGradDUZYZ = gradGradDUcompYZ.component(vector::Z);

volScalarField gradGradDUZZX = gradGradDUcompZX.component(vector::Z);
volScalarField gradGradDUZZY = gradGradDUcompZY.component(vector::Z);
volScalarField gradGradDUZZZ = gradGradDUcompZZ.component(vector::Z);

forAll(extrapGradDU.internalField(), facei)
{
    const label own = mesh.owner()[facei];
    const label nei = mesh.neighbour()[facei];
    const vector deltaOwn = mesh.Cf()[facei] - mesh.C()[own];
    const vector deltaNei = mesh.Cf()[facei] - mesh.C()[nei];
    const tensor& gradDUOwn = gradDU.internalField()[own];
    const tensor& gradDUNei = gradDU.internalField()[nei];

    // as there is there is no thirdOrderTensor class, we will
    // calculate (deltaOwn&gradGradDUOwn) out manually
    // tensor deltaOwnDotgradGradDUOwn = tensor::zero;
    // tensor deltaNeiDotgradGradDUNei = tensor::zero;

    // deltaOwnDotgradGradDUOwn[tensor::XX] = deltaOwn & gradGradDUcompXX.internalField()[own];
    // deltaNeiDotgradGradDUNei[tensor::XX] = deltaNei & gradGradDUcompXX.internalField()[nei];
    // deltaOwnDotgradGradDUOwn[tensor::XY] = deltaOwn & gradGradDUcompXY.internalField()[own];
    // deltaNeiDotgradGradDUNei[tensor::XY] = deltaNei & gradGradDUcompXY.internalField()[nei];
    // deltaOwnDotgradGradDUOwn[tensor::XZ] = deltaOwn & gradGradDUcompXZ.internalField()[own];
    // deltaNeiDotgradGradDUNei[tensor::XZ] = deltaNei & gradGradDUcompXZ.internalField()[nei];

    // deltaOwnDotgradGradDUOwn[tensor::YX] = deltaOwn & gradGradDUcompYX.internalField()[own];
    // deltaNeiDotgradGradDUNei[tensor::YX] = deltaNei & gradGradDUcompYX.internalField()[nei];
    // deltaOwnDotgradGradDUOwn[tensor::YY] = deltaOwn & gradGradDUcompYY.internalField()[own];
    // deltaNeiDotgradGradDUNei[tensor::YY] = deltaNei & gradGradDUcompYY.internalField()[nei];
    // deltaOwnDotgradGradDUOwn[tensor::YZ] = deltaOwn & gradGradDUcompYZ.internalField()[own];
    // deltaNeiDotgradGradDUNei[tensor::YZ] = deltaNei & gradGradDUcompYZ.internalField()[nei];

    // deltaOwnDotgradGradDUOwn[tensor::ZX] = deltaOwn & gradGradDUcompZX.internalField()[own];
    // deltaNeiDotgradGradDUNei[tensor::ZX] = deltaNei & gradGradDUcompZX.internalField()[nei];
    // deltaOwnDotgradGradDUOwn[tensor::ZY] = deltaOwn & gradGradDUcompZY.internalField()[own];
    // deltaNeiDotgradGradDUNei[tensor::ZY] = deltaNei & gradGradDUcompZY.internalField()[nei];
    // deltaOwnDotgradGradDUOwn[tensor::ZZ] = deltaOwn & gradGradDUcompZZ.internalField()[own];
    // deltaNeiDotgradGradDUNei[tensor::ZZ] = deltaNei & gradGradDUcompZZ.internalField()[nei];

    scalar gradGradDUXXXOwn = gradGradDUXXX.internalField()[own];
    scalar gradGradDUXXYOwn = gradGradDUXXY.internalField()[own];
    scalar gradGradDUXXZOwn = gradGradDUXXZ.internalField()[own];

    scalar gradGradDUXYXOwn = gradGradDUXYX.internalField()[own];
    scalar gradGradDUXYYOwn = gradGradDUXYY.internalField()[own];
    scalar gradGradDUXYZOwn = gradGradDUXYZ.internalField()[own];

    scalar gradGradDUXZXOwn = gradGradDUXZX.internalField()[own];
    scalar gradGradDUXZYOwn = gradGradDUXZY.internalField()[own];
    scalar gradGradDUXZZOwn = gradGradDUXZZ.internalField()[own];


    scalar gradGradDUYXXOwn = gradGradDUYXX.internalField()[own];
    scalar gradGradDUYXYOwn = gradGradDUYXY.internalField()[own];
    scalar gradGradDUYXZOwn = gradGradDUYXZ.internalField()[own];

    scalar gradGradDUYYXOwn = gradGradDUYYX.internalField()[own];
    scalar gradGradDUYYYOwn = gradGradDUYYY.internalField()[own];
    scalar gradGradDUYYZOwn = gradGradDUYYZ.internalField()[own];

    scalar gradGradDUYZXOwn = gradGradDUYZX.internalField()[own];
    scalar gradGradDUYZYOwn = gradGradDUYZY.internalField()[own];
    scalar gradGradDUYZZOwn = gradGradDUYZZ.internalField()[own];


    scalar gradGradDUZXXOwn = gradGradDUZXX.internalField()[own];
    scalar gradGradDUZXYOwn = gradGradDUZXY.internalField()[own];
    scalar gradGradDUZXZOwn = gradGradDUZXZ.internalField()[own];

    scalar gradGradDUZYXOwn = gradGradDUZYX.internalField()[own];
    scalar gradGradDUZYYOwn = gradGradDUZYY.internalField()[own];
    scalar gradGradDUZYZOwn = gradGradDUZYZ.internalField()[own];

    scalar gradGradDUZZXOwn = gradGradDUZZX.internalField()[own];
    scalar gradGradDUZZYOwn = gradGradDUZZY.internalField()[own];
    scalar gradGradDUZZZOwn = gradGradDUZZZ.internalField()[own];


    // nei
    scalar gradGradDUXXXNei = gradGradDUXXX.internalField()[nei];
    scalar gradGradDUXXYNei = gradGradDUXXY.internalField()[nei];
    scalar gradGradDUXXZNei = gradGradDUXXZ.internalField()[nei];

    scalar gradGradDUXYXNei = gradGradDUXYX.internalField()[nei];
    scalar gradGradDUXYYNei = gradGradDUXYY.internalField()[nei];
    scalar gradGradDUXYZNei = gradGradDUXYZ.internalField()[nei];

    scalar gradGradDUXZXNei = gradGradDUXZX.internalField()[nei];
    scalar gradGradDUXZYNei = gradGradDUXZY.internalField()[nei];
    scalar gradGradDUXZZNei = gradGradDUXZZ.internalField()[nei];


    scalar gradGradDUYXXNei = gradGradDUYXX.internalField()[nei];
    scalar gradGradDUYXYNei = gradGradDUYXY.internalField()[nei];
    scalar gradGradDUYXZNei = gradGradDUYXZ.internalField()[nei];

    scalar gradGradDUYYXNei = gradGradDUYYX.internalField()[nei];
    scalar gradGradDUYYYNei = gradGradDUYYY.internalField()[nei];
    scalar gradGradDUYYZNei = gradGradDUYYZ.internalField()[nei];

    scalar gradGradDUYZXNei = gradGradDUYZX.internalField()[nei];
    scalar gradGradDUYZYNei = gradGradDUYZY.internalField()[nei];
    scalar gradGradDUYZZNei = gradGradDUYZZ.internalField()[nei];


    scalar gradGradDUZXXNei = gradGradDUZXX.internalField()[nei];
    scalar gradGradDUZXYNei = gradGradDUZXY.internalField()[nei];
    scalar gradGradDUZXZNei = gradGradDUZXZ.internalField()[nei];

    scalar gradGradDUZYXNei = gradGradDUZYX.internalField()[nei];
    scalar gradGradDUZYYNei = gradGradDUZYY.internalField()[nei];
    scalar gradGradDUZYZNei = gradGradDUZYZ.internalField()[nei];

    scalar gradGradDUZZXNei = gradGradDUZZX.internalField()[nei];
    scalar gradGradDUZZYNei = gradGradDUZZY.internalField()[nei];
    scalar gradGradDUZZZNei = gradGradDUZZZ.internalField()[nei];

    tensor deltaOwnDotgradGradDUOwn = tensor
    (
        deltaOwn.x()*gradGradDUXXXOwn + deltaOwn.y()*gradGradDUYXXOwn + deltaOwn.z()*gradGradDUZXXOwn,
        deltaOwn.x()*gradGradDUXXYOwn + deltaOwn.y()*gradGradDUYXYOwn + deltaOwn.z()*gradGradDUZXYOwn,
        deltaOwn.x()*gradGradDUXXZOwn + deltaOwn.y()*gradGradDUYXZOwn + deltaOwn.z()*gradGradDUZXZOwn,

        deltaOwn.x()*gradGradDUXYXOwn + deltaOwn.y()*gradGradDUYYXOwn + deltaOwn.z()*gradGradDUZYXOwn,
        deltaOwn.x()*gradGradDUXYYOwn + deltaOwn.y()*gradGradDUYYYOwn + deltaOwn.z()*gradGradDUZYYOwn,
        deltaOwn.x()*gradGradDUXYZOwn + deltaOwn.y()*gradGradDUYYZOwn + deltaOwn.z()*gradGradDUZYZOwn,

        deltaOwn.x()*gradGradDUXZXOwn + deltaOwn.y()*gradGradDUYZXOwn + deltaOwn.z()*gradGradDUZZXOwn,
        deltaOwn.x()*gradGradDUXZYOwn + deltaOwn.y()*gradGradDUYZYOwn + deltaOwn.z()*gradGradDUZZYOwn,
        deltaOwn.x()*gradGradDUXZZOwn + deltaOwn.y()*gradGradDUYZZOwn + deltaOwn.z()*gradGradDUZZZOwn
    );

    tensor deltaNeiDotgradGradDUNei = tensor
    (
        deltaNei.x()*gradGradDUXXXNei + deltaNei.y()*gradGradDUYXXNei + deltaNei.z()*gradGradDUZXXNei,
        deltaNei.x()*gradGradDUXXYNei + deltaNei.y()*gradGradDUYXYNei + deltaNei.z()*gradGradDUZXYNei,
        deltaNei.x()*gradGradDUXXZNei + deltaNei.y()*gradGradDUYXZNei + deltaNei.z()*gradGradDUZXZNei,

        deltaNei.x()*gradGradDUXYXNei + deltaNei.y()*gradGradDUYYXNei + deltaNei.z()*gradGradDUZYXNei,
        deltaNei.x()*gradGradDUXYYNei + deltaNei.y()*gradGradDUYYYNei + deltaNei.z()*gradGradDUZYYNei,
        deltaNei.x()*gradGradDUXYZNei + deltaNei.y()*gradGradDUYYZNei + deltaNei.z()*gradGradDUZYZNei,

        deltaNei.x()*gradGradDUXZXNei + deltaNei.y()*gradGradDUYZXNei + deltaNei.z()*gradGradDUZZXNei,
        deltaNei.x()*gradGradDUXZYNei + deltaNei.y()*gradGradDUYZYNei + deltaNei.z()*gradGradDUZZYNei,
        deltaNei.x()*gradGradDUXZZNei + deltaNei.y()*gradGradDUYZZNei + deltaNei.z()*gradGradDUZZZNei
    );


    // get average of extrapolated values
    tensor extrapFaceGrad =
        0.5*
        (
            gradDUOwn + (deltaOwnDotgradGradDUOwn)
            +
            gradDUNei + (deltaNeiDotgradGradDUNei)
        );

    extrapGradDU.internalField()[facei] = extrapFaceGrad;

}

// correction is zero on boundary
forAll(extrapGradDU.boundaryField(), patchi)
{
    extrapGradDU.boundaryField()[patchi] =
        averageGradDU.boundaryField()[patchi];
}

// calculate thirdOrderTerm
volVectorField divThirdOrderTerm
(
    "thirdOrderTerm",
    fvc::div
    (
        (2*muf+lambdaf)*mesh.Sf() & (extrapGradDU - averageGradDU)
    )
);

// if(runTime.outputTime())
// {
//     divThirdOrderTerm.write();
//     averageGradDU.write();
//     extrapGradDU.write();
// }
