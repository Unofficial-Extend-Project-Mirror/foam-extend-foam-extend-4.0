if (!isType<zeroGradientFvPatchScalarField>(k_.boundaryField()[patchi]))
{
    FatalErrorIn("wall-function evaluation")
        << k_.boundaryField()[patchi].type()
        << " is the wrong k patchField type on patch "
        << curPatch.name() << nl
        << "    should be zeroGradient"
        << exit(FatalError);
}

if (!isType<zeroGradientFvPatchScalarField>(omega_.boundaryField()[patchi]))
{
    FatalErrorIn("wall-function evaluation")
        << omega_.boundaryField()[patchi].type()
        << " is the wrong omega patchField type on patch "
        << curPatch.name() << nl
        << "    should be zeroGradient"
        << exit(FatalError);
}
