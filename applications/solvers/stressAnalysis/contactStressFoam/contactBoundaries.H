{
    volVectorField::GeometricBoundaryField& Upatches = U.boundaryField();

    const volTensorField::GeometricBoundaryField& gradUpatches =
        gradU.boundaryField();

    const surfaceVectorField::GeometricBoundaryField& Apatches =
        mesh.Sf().boundaryField();
    const surfaceScalarField::GeometricBoundaryField& magApatches =
        mesh.magSf().boundaryField();

    vectorField nGradPatch = Apatches[gradPatch]/magApatches[gradPatch];

    vectorField nDirPatch = Apatches[dirPatch]/magApatches[dirPatch];

    // Update contact
    contactPair.updateContact(U);
    reversePair.updateContact(U);

    const scalarField& touchFraction = contactPair.touchFraction();

    const scalarField& reverseFraction = reversePair.touchFraction();

//Info << "touchFraction: " << touchFraction << endl;

    // Mark contact surfaces
    contactArea.boundaryField()[dirPatch] = touchFraction;
    contactArea.boundaryField()[gradPatch] = reverseFraction;

    // Cast will fall over in incorrect
    directionMixedFvPatchVectorField& UpatchDir =
        refCast<directionMixedFvPatchVectorField>(Upatches[dirPatch]);

    // set the traction and value for directionMixed patch
    UpatchDir.valueFraction() =
        (1.0 - urf)*UpatchDir.valueFraction()
      + I*urf*touchFraction;

//     UpatchDir.value() = contactPair.slaveDisplacement();

    UpatchDir.refValue() =
        nDirPatch*
        min
        (
            (nDirPatch & UpatchDir) + urf*touchFraction*touchTolerance,
            (nDirPatch & contactPair.slaveDisplacement())
        );

    // traction[dirPatch] = 0 because there's no friction!!

    // set the traction for the gradient patch (using the new valueFraction!)

    vectorField newTraction =
        nGradPatch*
        (
            reversePair.slavePressure
            (
                nDirPatch &
                (
                    mu.value()*
                    (
                        gradUpatches[dirPatch]
                      + gradUpatches[dirPatch].T()
                    )
                  + I*(lambda.value()*tr(gradUpatches[dirPatch]))
                ) & nDirPatch
            )
        );

//     Info<< "contact: "
//         <<
//             (nDirPatch &
//                 (
//                     mu.value()*
//                     (
//                         gradUpatches[dirPatch]
//                       + gradUpatches[dirPatch].T()
//                     )
//                   + I*(lambda.value()*tr(gradUpatches[dirPatch]))
//                 ) & nDirPatch)*rho.value()
//         << endl;

    traction[gradPatch] =
        (1.0 - urf)*traction[gradPatch]
      + urf*newTraction*reverseFraction;
}
