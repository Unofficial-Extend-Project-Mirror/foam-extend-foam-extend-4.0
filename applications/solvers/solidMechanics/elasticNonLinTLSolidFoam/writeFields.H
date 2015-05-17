if (runTime.outputTime())
{
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

    Info<< "Max epsilonEq = " << max(epsilonEq).value()
        << endl;

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

    Info<< "Max sigmaEq = " << max(sigmaEq).value()
        << endl;

    //- Calculate Cauchy stress
    volTensorField F = I + gradU;
    volScalarField J = det(F);

    //- update density
    rho = rho/J;

    volSymmTensorField sigmaCauchy
    (
        IOobject
        (
           "sigmaCauchy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (1/J) * symm(F.T() & sigma & F)
    );

    //- Cauchy von Mises stress
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
        sqrt((3.0/2.0)*magSqr(dev(sigmaCauchy)))
    );

    Info<< "Max sigmaCauchyEq = " << max(sigmaCauchyEq).value()
        << endl;


    volTensorField Finv = inv(F);
    volSymmTensorField epsilonAlmansi
    (
        IOobject
        (
            "epsilonAlmansi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        symm(Finv & epsilon & Finv.T())
    );

//     volVectorField traction
//       (
//        IOobject
//        (
//     "traction",
//     runTime.timeName(),
//     mesh,
//     IOobject::NO_READ,
//     IOobject::AUTO_WRITE
//     ),
//        mesh,
//        dimensionedVector("zero", dimForce/dimArea, vector::zero),
//        calculatedFvPatchVectorField::typeName
//        );
//     forAll(traction.boundaryField(), patchi)
//       {
//     const tensorField& Fbinv = Finv.boundaryField()[patchi];
//     vectorField nCurrent = Fbinv & n.boundaryField()[patchi];
//     traction.boundaryField()[patchi] =
//       nCurrent & sigmaCauchy.boundaryField()[patchi];
//       }

//     //- write boundary forces
//     //- integrate (sigma2PK & F) over reference area
//     //- which is equivalent to integrating sigmaCauchy
//     //- over the deformed area
//     Info << nl;
//     forAll(mesh.boundary(), patchi)
//       {
//     Info << "Patch " << mesh.boundary()[patchi].name() << endl;
//     const tensorField& Fb = F.boundaryField()[patchi];
//     vectorField totalForce = mesh.Sf().boundaryField()[patchi] & (sigma.boundaryField()[patchi] & Fb);
//     //vectorField totalForce2 = Sf.boundaryField()[patchi] & (sigmaCauchy.boundaryField()[patchi]);

//     vector force = sum( totalForce );
//     //vector force2 = sum( totalForce2 );
//     Info << "\ttotal force is " << force << " N" << endl;
//     //Info << "\ttotal force2 is " << force2 << " N" << endl;

//     const tensorField& Fbinv = Finv.boundaryField()[patchi];
//     vectorField nCurrent = Fbinv & n.boundaryField()[patchi];
//     nCurrent /= mag(nCurrent);
//     scalar normalForce = sum( nCurrent & totalForce );
//     Info << "\tnormal force is " << normalForce << " N" << endl;
//     scalar shearForce = mag(sum( (I - sqr(nCurrent)) & totalForce ));
//     Info << "\tshear force is " << shearForce << " N" << endl;

    //if(mesh.boundary()[patchi].name() == "right")
    //{
    //const vectorField& nOrig = n.boundaryField()[patchi];
    //Info << "\tNormal force on right is " << (nCurrent & totalForce) << nl << endl;
    //Info << "\tShear force on right is " << ((I - sqr(nCurrent)) & totalForce) << nl << endl;
    //Info << "\tpatch gradient is " << U.boundaryField()[patchi].snGrad() << endl;
    //Info << "\tpatch gradient (norm) is " << (nCurrent & U.boundaryField()[patchi].snGrad()) << endl;
    //Info << "\tpatch gradient (shear) is " << ((I - sqr(nCurrent)) & U.boundaryField()[patchi].snGrad()) << endl;
    //Info << "\tpatch Almansi (normal) is " << (nCurrent & (nCurrent & epsilonAlmansi.boundaryField()[patchi])) << endl;
    //Info << "\tpatch Almansi (shear) is " << ( (I - sqr(nCurrent)) & (nCurrent & epsilonAlmansi.boundaryField()[patchi])) << endl;
    //Info << "\tpatch Green (normal) is " << (nOrig & (nOrig & epsilon.boundaryField()[patchi])) << endl;
    //Info << "\tpatch Green (shear) is " << ( (I - sqr(nOrig)) & (nOrig & epsilon.boundaryField()[patchi])) << endl;
    //Info << "\tpatch Cauchy stress (normal) is " << (nCurrent & (nCurrent & sigmaCauchy.boundaryField()[patchi])) << endl;
    //}

//     if(mesh.boundary()[patchi].type() != "empty")
//       {
//             vector Sf0 = Sf.boundaryField()[patchi][0];
//             symmTensor sigma0 = sigmaCauchy.boundaryField()[patchi][0];
//             Info << "sigmab[0] is " << sigma0 << nl
//                  << "Sfb is " << Sf0 << nl
//                  << "force is " << (Sf.boundaryField()[patchi][0]&sigma.boundaryField()[patchi][0]) << nl
//          << "Sfx*sigmaxx " << (Sf0[vector::X]*sigma0[symmTensor::XX]) <<nl
//                  << "Sfy*sigmaxy " << (Sf0[vector::Y]*sigma0[symmTensor::XY]) << nl
//                  << "Sfx*sigmayx " << (Sf0[vector::X]*sigma0[symmTensor::XY]) << nl
//                  << "Sfy*sigmayy " << (Sf0[vector::Y]*sigma0[symmTensor::YY]) << nl
//                  << endl;
//       }
//    Info << endl;
//    }

    runTime.write();
}
