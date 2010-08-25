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

volScalarField scalarField
(
    fieldObject,
    *secondMeshPtr
);

const cellShapeList& cellShapes = secondMeshPtr->cellShapes();

// hexa's
if (which_type == Z_HEX08)
{
    const cellModel& hex = *(cellModeller::lookup("hex"));
    //const cellModel& wedge = *(cellModeller::lookup("wedge"));

    label counter = 1;
    for (label celli=0; celli<nCells; celli++)
    {
        const cellShape& cellShape = cellShapes[celli];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == hex) // || (cellModel == wedge))
        {
            var_array[counter++] = scalarField[celli];
        }
    }
}

// penta's
if (which_type == Z_PEN06)
{
    const cellModel& prism = *(cellModeller::lookup("prism"));

    label counter = 1;
    for (label n=0; n<nCells; n++)
    {
        const cellShape& cellShape = cellShapes[n];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == prism)
        {
            var_array[counter++] = scalarField[n];
        }
    }
}

// pyramids's
if (which_type == Z_PYR05)
{
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
            
    label counter = 1;
    for (label n=0; n<nCells; n++)
    {
        const cellShape& cellShape = cellShapes[n];
        const cellModel& cellModel = cellShape.model();
        
        if (cellModel == pyr)
        {
            var_array[counter++] = scalarField[n];
        }
    }
}

// tet's
if (which_type == Z_TET04)
{
    const cellModel& tet = *(cellModeller::lookup("tet"));

    label counter = 1;
    for (label n=0; n<nCells; n++)
    {
        const cellShape& cellShape = cellShapes[n];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == tet)
        {
            var_array[counter++] = scalarField[n];
        }
    }
}

if (which_type == Z_NFACED)
{
    const cellList& cells = secondMeshPtr->cells();

    label counter = 1;
    for (label n=0; n<nCells; n++)
    {
        const labelList& points = cellShapes[n];
        label nFacesInCell = cells[n].size();

        if ((nFacesInCell == 6) && (points.size() == 8))
        {}
        else if ((nFacesInCell == 4) && (points.size() == 4))
        {}
        else if (nFacesInCell == 5)
        {
            if (points.size() == 6)
            {}
            else if (points.size() == 5)
            {}
            else
            {
                var_array[counter++] = scalarField[n];
            }
        }
        else
        {
            var_array[counter++] = scalarField[n];
        }
    }

}
