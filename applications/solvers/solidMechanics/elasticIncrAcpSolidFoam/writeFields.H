
if (runTime.outputTime() || topoChange)
{
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

//     Info << "\nCalculate maximal principal stress ..." << flush;
//     // Principal stresses
//     volVectorField sigmaMax
//     (
//         IOobject
//         (
//             "sigmaMax",
//             runTime.timeName(),
//             mesh,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         mesh,
//         dimensionedVector("sigmaMax", dimPressure, vector::zero)
//     );
//     vectorField& sigmaMaxI = sigmaMax.internalField();

//     forAll (sigmaMaxI, cellI)
//     {
//         vector eValues = eigenValues(sigma.internalField()[cellI]);
//         tensor eVectors = eigenVectors(sigma.internalField()[cellI]);

//         scalar maxEValue = 0;
//         label iMax = -1;
//         forAll(eValues, i)
//         {
//             if (eValues[i] > maxEValue)
//             {
//                 maxEValue = eValues[i];
//                 iMax = i;
//             }
//         }

//         if (iMax != -1)
//         {
//             if (iMax == 0)
//             {
//                 sigmaMaxI[cellI] = eVectors.x()*eValues.x();
//             }
//             else if (iMax == 1)
//             {
//                 sigmaMaxI[cellI] = eVectors.y()*eValues.y();
//             }
//             else if (iMax == 2)
//             {
//                 sigmaMaxI[cellI] = eVectors.z()*eValues.z();
//             }
//         }
//     }

//     //- boundary traction
//     volVectorField tractionBoundary
//     (
//         IOobject
//         (
//             "tractionBoundary",
//             runTime.timeName(),
//             mesh,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         mesh,
//         dimensionedVector("zero", dimForce/dimArea, vector::zero)
//     );
//     surfaceVectorField n = mesh.Sf()/mesh.magSf();
//     forAll(tractionBoundary.boundaryField(), patchi)
//     {
//     if(mesh.boundaryMesh()[patchi].type() != processorPolyPatch::typeName)
//     {
//         tractionBoundary.boundaryField()[patchi] =
//         n.boundaryField()[patchi] & sigma.boundaryField()[patchi];
//     }
// }


    //- cohesive damage and cracking, and GII and GII
    volScalarField damageAndCracks
    (
        IOobject
        (
            "damageAndCracks",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0),
        calculatedFvPatchVectorField::typeName
    );

    volScalarField GI
    (
        IOobject
        (
            "GI",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       mesh,
       dimensionedScalar("zero", dimless, 0.0),
       calculatedFvPatchVectorField::typeName
    );

    volScalarField GII
    (
        IOobject
        (
             "GII",
             runTime.timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0),
        calculatedFvPatchVectorField::typeName
    );

    forAll(DU.boundaryField(), patchi)
    {
        // if(DU.boundaryField()[patchi].type() == cohesiveLawMultiMatFvPatchVectorField::typeName)
        if(DU.boundaryField()[patchi].type() == solidCohesiveFvPatchVectorField::typeName)
        {
            // cohesiveLawMultiMatFvPatchVectorField& DUpatch =
            //   refCast<cohesiveLawMultiMatFvPatchVectorField>(DU.boundaryField()[patchi]);
            solidCohesiveFvPatchVectorField& DUpatch =
                refCast<solidCohesiveFvPatchVectorField>(DU.boundaryField()[patchi]);

            GI.boundaryField()[patchi] = DUpatch.GI();
            GII.boundaryField()[patchi] = DUpatch.GII();
            damageAndCracks.boundaryField()[patchi] = DUpatch.crackingAndDamage();
        }
    }

    //Info << "done" << endl;

    runTime.writeNow();
}
