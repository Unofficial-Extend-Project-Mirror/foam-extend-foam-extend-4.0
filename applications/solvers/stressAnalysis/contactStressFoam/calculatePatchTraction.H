    Info << "Calculating patch traction" << endl;
    PtrList<vectorField> traction(mesh.boundary().size());

    {
        const surfaceVectorField::GeometricBoundaryField& Apatches =
            mesh.Sf().boundaryField();

        const surfaceScalarField::GeometricBoundaryField& magApatches =
            mesh.magSf().boundaryField();

        forAll(traction, patchI)
        {
            vectorField n = Apatches[patchI]/magApatches[patchI];

            traction.set
            (
                patchI,
                new vectorField
                (
                    patchTraction[patchI]
                  - patchPressure[patchI]*n
                )
            );
        }
    }

