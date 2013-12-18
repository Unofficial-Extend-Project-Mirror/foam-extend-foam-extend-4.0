    if (runTime.outputTime())
    {
        volScalarField sigmaCauchyEq
        (
            IOobject
            (
                "sigmaCauchyEq",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqrt((3.0/2.0)*magSqr(dev(sigma)))
        );

        Info<< "Max sigmaCauchyEq = " << max(sigmaCauchyEq).value()
            << endl;

//         volScalarField sigmaCauchyxx
//         (
//             IOobject
//             (
//                 "sigmaCauchyxx",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             sigma.component(symmTensor::XX)
//         );

//         Info<< "Max sigmaCauchyxx = " << max(sigmaCauchyxx).value()
//             << endl;

//         volScalarField sigmaCauchyyy
//         (
//             IOobject
//             (
//                 "sigmaCauchyyy",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             sigma.component(symmTensor::YY)
//         );

//         Info<< "Max sigmaCauchyyy = " << max(sigmaCauchyyy).value()
//             << endl;

//         volScalarField sigmaCauchyzz
//         (
//             IOobject
//             (
//                 "sigmaCauchyzz",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             sigma.component(symmTensor::ZZ)
//         );

//         Info<< "Max sigmaCauchyzz = " << max(sigmaCauchyzz).value()
//             << endl;

//         volScalarField sigmaCauchyxy
//         (
//             IOobject
//             (
//                 "sigmaCauchyxy",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             sigma.component(symmTensor::XY)
//         );

//         volScalarField sigmaCauchyxz
//         (
//             IOobject
//             (
//                 "sigmaCauchyxz",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             sigma.component(symmTensor::XZ)
//         );

//         volScalarField sigmaCauchyyz
//         (
//             IOobject
//             (
//                 "sigmaCauchyyz",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             sigma.component(symmTensor::YZ)
//         );

        runTime.write();
    }
