{
    volVectorField::GeometricBoundaryField& Upatches = U.boundaryField();

    const volTensorField::GeometricBoundaryField& gradUpatches =
        gradU.boundaryField();

    const surfaceVectorField::GeometricBoundaryField& Apatches =
        mesh.Sf().boundaryField();

    const surfaceScalarField::GeometricBoundaryField& magApatches =
        mesh.magSf().boundaryField();

    forAll(Upatches, patchI)
    {
        if
        (
            Upatches[patchI].type()
         == fixedGradientFvPatchVectorField::typeName
        )
        {
            fixedGradientFvPatchVectorField& Upatch =
                refCast<fixedGradientFvPatchVectorField>(Upatches[patchI]);

            vectorField n = Apatches[patchI]/magApatches[patchI];

            Upatch.gradient() =
#           include "traction.H"
        }

        if
        (
            Upatches[patchI].type()
         == directionMixedFvPatchVectorField::typeName
        )
        {

            directionMixedFvPatchVectorField& Upatch =
                refCast<directionMixedFvPatchVectorField>(Upatches[patchI]);

            vectorField n = Apatches[patchI]/magApatches[patchI];

            Upatch.refGrad() =
#           include "traction.H"
       }
    }
}
