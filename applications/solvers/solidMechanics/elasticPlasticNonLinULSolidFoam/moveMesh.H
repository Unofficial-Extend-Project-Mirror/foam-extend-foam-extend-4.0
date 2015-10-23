if(moveMeshMethod == "inverseDistance")
{
#   include "moveMeshInverseDistance.H"
}
else if(moveMeshMethod == "leastSquares")
{
#   include "moveMeshLeastSquares.H"
}
else
{
    FatalErrorIn(args.executable())
        << "move mesh method " << moveMeshMethod << " not recognised" << nl
        << "available methods are:" << nl
        << "inverseDistance" << nl
        << "leastSquares" << exit(FatalError);
}
