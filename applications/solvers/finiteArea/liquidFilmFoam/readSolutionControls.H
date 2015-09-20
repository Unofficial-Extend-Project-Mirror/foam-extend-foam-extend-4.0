int nCorr = 0;

if (aMesh.solutionDict().found("nOuterCorrectors"))
{
    nCorr =
        readInt(aMesh.solutionDict().lookup("nOuterCorrectors"));
}
else
{
    FatalErrorIn(args.executable())
        << "Cannot find number of correctors"
        << abort(FatalError);
}
