for (int n1=startTime; n1<endTime; n1++)
{
    if (fieldName.size() > 2 && fieldName(fieldName.size() - 2, 2) == "_0")
    {
        variableGood = false;
        break;
    }

    IOobject fieldObjectHeader
    (
        fieldName,
        Times[n1].name(),
        mesh,
        IOobject::NO_READ
    );

    variableGood = variableGood && fieldObjectHeader.headerOk();
}
