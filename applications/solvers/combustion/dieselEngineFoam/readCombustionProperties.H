Info<< "Reading combustion properties\n" << endl;

IOdictionary combustionProperties
(
    IOobject
    (
        "combustionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar Cmix
(
    combustionProperties.lookup("Cmix")
);
