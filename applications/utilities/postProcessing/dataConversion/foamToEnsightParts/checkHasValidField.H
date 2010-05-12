// check that the variable is present for all times
//
bool hasValidField = true;
{
    for (label i=0; i < timeDirs.size() && hasValidField; ++i)
    {
        if (fieldName.size() > 2 && fieldName(fieldName.size() - 2, 2) == "_0")
        {
            hasValidField = false;
            break;
        }

        IOobject ioHeader
        (
            fieldName,
            timeDirs[i].name(),
            mesh,
            IOobject::NO_READ
        );

        hasValidField = ioHeader.headerOk();
    }
}
