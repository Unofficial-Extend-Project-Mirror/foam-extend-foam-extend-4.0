{
    forAll(mesh.boundary(), patchID)
    {
        if(U.boundaryField()[patchID].type()
     == solidDirectionMixedFvPatchVectorField::typeName
    )
    {
        solidDirectionMixedFvPatchVectorField& loadingPatch =
            refCast<solidDirectionMixedFvPatchVectorField>
            (
                U.boundaryField()[patchID]
            );

            tensorField Finv = inv(I + gradU);
            vectorField newN = Finv & n.boundaryField()[patchID];
            newN /= mag(newN);
            loadingPatch.valueFraction() = sqr(newN);

            //- set gradient
            loadingPatch.refGrad() =
            (
                //Traction
                ( (mu.boundaryField()[patchID] + lambda.boundaryField()[patchID]) * (n.boundaryField()[patchID] & gradU.boundaryField()[patchID]) )
              - ( mu.boundaryField()[patchID] * (n.boundaryField()[patchID] & gradU.boundaryField()[patchID].T()) )
              - ( mu.boundaryField()[patchID] * ( n.boundaryField()[patchID] & (gradU.boundaryField()[patchID] & gradU.boundaryField()[patchID].T()) ) )
              - ( lambda.boundaryField()[patchID] * tr(gradU.boundaryField()[patchID]) * n.boundaryField()[patchID] )
              - ( 0.5 * lambda.boundaryField()[patchID] * tr(gradU.boundaryField()[patchID] & gradU.boundaryField()[patchID].T()) * n.boundaryField()[patchID] )
            )
            /
            (2.0*mu.boundaryField()[patchID] + lambda.boundaryField()[patchID]);
        }
    }
}
