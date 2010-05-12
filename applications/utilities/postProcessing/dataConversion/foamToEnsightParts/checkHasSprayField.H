// check that the spray variable is present for this time
//
bool hasSprayField = true;
{
    IOobject ioHeader
    (
        fieldName,
        mesh.time().timeName(),
        "lagrangian",
        mesh,
        IOobject::NO_READ
    );

    hasSprayField = ioHeader.headerOk();
}
