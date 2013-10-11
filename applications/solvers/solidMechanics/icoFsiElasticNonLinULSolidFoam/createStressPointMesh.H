
// Create point stress mesh for interpolation
pointMesh pStressMesh(stressMesh);

pointPatchInterpolation patchPointInterpolator(stressMesh);


// Create point displacement field (check fixedValue patches)

wordList types
(
    pStressMesh.boundary().size(), 
    calculatedFvPatchVectorField::typeName
);

// wordList types = DU.boundaryField().types();

forAll(DU.boundaryField().types(), patchI)
{
    if 
    (
        DU.boundaryField().types()[patchI] 
     == fixedValueFvPatchVectorField::typeName
    )
    {
        types[patchI] = fixedValueFvPatchVectorField::typeName;
    }
}

pointVectorField pointDU
(
    IOobject
    (
        "pointDU",
        runTime.timeName(),
        stressMesh
    ),
    pStressMesh,
    dimensionedVector("zero", dimLength, vector::zero),
    types
);
