vector latticePosition;

vector globalPosition;

if (latticeStructure == "SC")
{
#   include "SC.H"
}

else if (latticeStructure == "FCC")
{
#   include "FCC.H"
}

else if (latticeStructure == "BCC")
{
#   include "BCC.H"
}

else
{
    FatalErrorIn("createPositions.H\n")
        << "latticeStructure " << latticeStructure
        << " not supported."
        << abort(FatalError);
}
