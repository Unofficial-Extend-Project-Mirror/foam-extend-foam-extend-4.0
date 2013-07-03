Info << "s1" << endl;

    fvMesh solidMesh
    (
        IOobject
        (
            "solid",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

Info << "s2" << endl;

