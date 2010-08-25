label offset = Num_variables - nFaMeshVariables;

if (nVar < offset || nVar >= Num_variables)
{
    return Z_UNDEF;
}

IOobject fieldObjectPtr
(
    faMeshFieldNames[faMeshVar2field[nVar-offset]],
    runTime.timeName(),
    *meshPtr,
    IOobject::NO_READ
);

if (!fieldObjectPtr.headerOk())
{
    return Z_UNDEF;
}

IOobject fieldObject
(
    faMeshFieldNames[faMeshVar2field[nVar-offset]],
    runTime.timeName(),
    *meshPtr,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
);

areaScalarField scalarField
(
    fieldObject,
    *faMeshPtr
);


const faceList& faces = faMeshPtr->faces();

label nFaces = faces.size();

if (which_type == Z_TRI03)
{
    label counter = 1;
    for (label facei=0; facei<nFaces; facei++)
    {
        label nPoints = faces[facei].size();
        if (nPoints == 3)
        {
            var_array[counter++] = scalarField[facei];
        }
    }
}

if (which_type == Z_QUA04)
{
    label counter = 1;
    for (label facei=0; facei<nFaces; facei++)
    {
        label nPoints = faces[facei].size();
        if (nPoints == 4)
        {
            var_array[counter++] = scalarField[facei];
        }
    }
}

if (which_type == Z_NSIDED)
{
    label counter = 1;
    for (label facei=0; facei<nFaces; facei++)
    {
        label nPoints = faces[facei].size();
        if ((nPoints != 3) && (nPoints != 4))
        {
            var_array[counter++] = scalarField[facei];
        }
    }
}
