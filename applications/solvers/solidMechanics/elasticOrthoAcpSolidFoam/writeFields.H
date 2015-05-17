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
    forAll(U.boundaryField(), patchi)
    {
        // if(U.boundaryField()[patchi].type() == cohesiveLawMultiMatFvPatchVectorField::typeName)
        if(U.boundaryField()[patchi].type() == solidCohesiveFvPatchVectorField::typeName)
        {
            // cohesiveLawMultiMatFvPatchVectorField& Upatch =
            //   refCast<cohesiveLawMultiMatFvPatchVectorField>(U.boundaryField()[patchi]);
            solidCohesiveFvPatchVectorField& Upatch =
                refCast<solidCohesiveFvPatchVectorField>(U.boundaryField()[patchi]);

            GI.boundaryField()[patchi] = Upatch.GI();
            GII.boundaryField()[patchi] = Upatch.GII();
            damageAndCracks.boundaryField()[patchi] = Upatch.crackingAndDamage();
        }
    }
    volScalarField GTotal("GTotal", GI + GII);
    GTotal.write();


    runTime.writeNow();
}
