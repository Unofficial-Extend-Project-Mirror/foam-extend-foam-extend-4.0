// check for lagrangian/positions information in the final directory

bool hasLagrangian = false;
if (timeDirs.size() > 1)
{
    IOobject io
    (
        "positions",
        timeDirs[timeDirs.size() - 1].name(),
        "lagrangian",
        mesh,
        IOobject::NO_READ
    );

    if (io.headerOk())
    {
        hasLagrangian = true;
    }
}
