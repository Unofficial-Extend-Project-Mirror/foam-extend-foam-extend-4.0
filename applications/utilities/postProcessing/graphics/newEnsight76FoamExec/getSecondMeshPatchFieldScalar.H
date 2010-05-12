
label secondMeshPartNum = nPatches+2;

if(sprayPtr)
{
    secondMeshPartNum += 1;
}

label patchi = which_part - secondMeshPartNum - 1;

label offset = Num_variables - nSecondMeshVariables - nFaMeshVariables;

if (nVar < offset || nVar >= (Num_variables-nFaMeshVariables))
{
    return Z_UNDEF;
}

IOobject fieldObjectPtr
(
    secondMeshFieldNames[secondMeshVar2field[nVar-offset]],
    runTime.timeName(),
    *secondMeshPtr,
    IOobject::NO_READ
);

if (!fieldObjectPtr.headerOk())
{
    return Z_UNDEF;
}

IOobject fieldObject
(
    secondMeshFieldNames[secondMeshVar2field[nVar-offset]],
    runTime.timeName(),
    *secondMeshPtr,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
);

volScalarField sf
(
    fieldObject,
    *secondMeshPtr
);

const scalarField& sfb = sf.boundaryField()[patchi];
const polyBoundaryMesh& bMesh = secondMeshPtr->boundaryMesh();

if (which_type == Z_TRI03)
{
    label counter = 1;
    for (label facei=0; facei<sfb.size(); facei++)
    {
        label nPoints = bMesh[patchi][facei].size();
        if (nPoints == 3)
        {
            var_array[counter++] = sfb[facei];
        }
    }
}

if (which_type == Z_QUA04)
{
    label counter = 1;
    for (label facei=0; facei<sfb.size(); facei++)
    {
        label nPoints = bMesh[patchi][facei].size();
        if (nPoints == 4)
        {
            var_array[counter++] = sfb[facei];
        }
    }
}

if (which_type == Z_NSIDED)
{
    label counter = 1;
    for (label facei=0; facei<sfb.size(); facei++)
    {
        label nPoints = bMesh[patchi][facei].size();
        if ((nPoints != 3) && (nPoints != 4))
        {
            var_array[counter++] = sfb[facei];
        }
    }
}
