forAll(p_.boundaryField(), patchI)
{
    if
    (
        p_.boundaryField()[patchI].type()
     == fixedGradientFvPatchScalarField::typeName
    )
    {
        Pout << "Seting interfac pressure gradient" << endl;

        fixedGradientFvPatchScalarField& pInterface =
            refCast<fixedGradientFvPatchScalarField>
            (
                p_.boundaryField()[patchI]
            );

        const vectorField n = mesh.boundary()[patchI].nf();

        pInterface.gradient() =
            (n&U_.boundaryField()[patchI])
           /rUA.boundaryField()[patchI]
          - phi_.boundaryField()[patchI]
           /mesh.magSf().boundaryField()[patchI]
           /rUA.boundaryField()[patchI]
          + nGradUn_
           /mesh.deltaCoeffs().boundaryField()[patchI]
           /rUA.boundaryField()[patchI];

        phi_.boundaryField()[patchI] +=
            pInterface.gradient()
           *rUA.boundaryField()[patchI]
           *mesh.magSf().boundaryField()[patchI];
    }
}
