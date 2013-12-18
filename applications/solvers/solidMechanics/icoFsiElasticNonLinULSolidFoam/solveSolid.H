
if (solidDdtScheme == fv::EulerDdtScheme<vector>::typeName)
{
#   include "solveSolidEuler.H"
}
else if (solidDdtScheme == fv::backwardDdtScheme<vector>::typeName)
{
#   include "solveSolidBackward.H"
}
else
{
    FatalErrorIn(args.executable())
        << "Wrong temporal (ddt) scheme for solid solver. "
            << "Possible schemes are: "
            << fv::EulerDdtScheme<vector>::typeName << " and "
            << fv::backwardDdtScheme<vector>::typeName
            << abort(FatalError);
}

