forAll(p().boundaryField(), patchI)
{
    if
    (
        (
            isA<elasticWallPressureFvPatchScalarField>(p().boundaryField()[patchI])
         && isA<elasticSlipWallVelocityFvPatchVectorField>(U().boundaryField()[patchI])
        )
     || (
            isA<elasticWallPressureFvPatchScalarField>(p().boundaryField()[patchI])
         && isA<elasticWallVelocityFvPatchVectorField>(U().boundaryField()[patchI])
        )
    )
    {
        word ddtScheme
        (
            mesh.schemesDict().ddtScheme("ddt(" + U().name() +')')
        );

        if
        (
            ddtScheme
         == fv::EulerDdtScheme<vector>::typeName
        )
        {
            phi().boundaryField()[patchI] = phi().oldTime().boundaryField()[patchI];
            rAUf.boundaryField()[patchI] = runTime().deltaT().value();
        }
        else if
        (
            ddtScheme
         == fv::backwardDdtScheme<vector>::typeName
        )
        {
            if(runTime().timeIndex() == 1)
            {
                phi().boundaryField()[patchI] = phi().oldTime().boundaryField()[patchI];
                rAUf.boundaryField()[patchI] = runTime().deltaT().value();

                phi().oldTime().oldTime();
            }
            else
            {
                scalar deltaT = runTime().deltaT().value();
                scalar deltaT0 = runTime().deltaT0().value();

                scalar Cn = 1 + deltaT/(deltaT + deltaT0);
                scalar Coo = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
                scalar Co = Cn + Coo;

                phi().boundaryField()[patchI] =
                    (Co/Cn)*phi().oldTime().boundaryField()[patchI]
                  - (Coo/Cn)*phi().oldTime().oldTime().boundaryField()[patchI];

                rAUf.boundaryField()[patchI] = runTime().deltaT().value()/Cn;
            }
        }
    }
}
