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

    //- boundary traction
    // volVectorField traction
    //   (
    //    IOobject
    //    (
    //     "traction",
    //     runTime.timeName(),
    //     mesh,
    //     IOobject::NO_READ,
    //     IOobject::AUTO_WRITE
    //     ),
    //    mesh,
    //    dimensionedVector("zero", dimForce/dimArea, vector::zero)
    //    );
    // forAll(traction.boundaryField(), patchi)
    //   {
    //     traction.boundaryField()[patchi] =
    //       n.boundaryField()[patchi] & sigma.boundaryField()[patchi];
    //   }

    // contact fields
    // if(contactBoundaries)
    //   {
    //     volScalarField stickSlipFaces
    //       (
    //        IOobject
    //        (
    //         "stickSlipFaces",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //         ),
    //        mesh,
    //        dimensionedScalar("zero", dimless, 0.0)
    //        );
    //     forAll(stickSlipFaces.boundaryField(), patchi)
    //       {
    //         if(U.boundaryField()[patchi].type()
    //            ==
    //            solidContactFvPatchVectorField::typeName)
    //           {
    //         const solidContactFvPatchVectorField& Upatch =
    //           refCast<const solidContactFvPatchVectorField>
    //           (U.boundaryField()[patchi]);

    //         if(!Upatch.master())
    //           {
    //             stickSlipFaces.boundaryField()[patchi] =
    //               Upatch.frictionContactModelPtr()->stickSlipFaces();
    //           }
    //           }
    //       }
    //     stickSlipFaces.write();
    //   }


    //- boundary forces
//     Info << nl;
//     forAll(mesh.boundary(), patchi)
//       {
//     Info << "Patch " << mesh.boundary()[patchi].name() << endl;
//     vectorField totalForce = mesh.Sf().boundaryField()[patchi] & sigma.boundaryField()[patchi];

//     vector force = sum( totalForce );
//     Info << "\ttotal force is " << force << " N" << endl;

//     const vectorField& nb = n.boundaryField()[patchi];
//     scalar normalForce = sum( nb & totalForce );
//     Info << "\tnormal force is " << normalForce << " N" << endl;
//     scalar shearForce = mag(sum( (I - sqr(nb)) & totalForce ));
//     Info << "\tshear force is " << shearForce << " N" << endl;
//     }

    runTime.write();
  }
