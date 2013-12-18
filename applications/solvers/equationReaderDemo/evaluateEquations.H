    // Update dependents
    eqnSrcR == turbulence->R();

    // scalars
    scalarOut = eqns.evaluateScalar("sOut");

    dimensionedScalarOut = eqns.evaluateDimensionedScalar("dsOut");

    eqns.evaluateDimensionedScalarField(dimensionedScalarFieldOut, "dsfOut");

    eqns.evaluateGeometricScalarField(volScalarFieldOut, "volSfOut");

    // vectors
    vectorOut.x() = eqns.evaluateScalar("vOut.x");
    vectorOut.y() = eqns.evaluateScalar("vOut.y");
    vectorOut.z() = eqns.evaluateScalar("vOut.z");

    dimensionedVectorOut.value().x() = eqns.evaluateScalar("dvOut.x");
    dimensionedVectorOut.value().y() = eqns.evaluateScalar("dvOut.y");
    dimensionedVectorOut.value().z() = eqns.evaluateScalar("dvOut.z");
    // There is currently no elegant way to check the dimensions of each
    // component of a "dimensionedType".  (dimensionedTypeFields and
    // GeometricTypeFields work automatically, though.)  This is how we do
    // it:
    dimensionedVectorOut.dimensions() = eqns.evaluateDimensions("dvOut.x");
    dimensionedVectorOut.dimensions() = eqns.evaluateDimensions("dvOut.y");
    dimensionedVectorOut.dimensions() = eqns.evaluateDimensions("dvOut.z");

    eqns.evaluateDimensionedTypeField
    (
        dimensionedVectorFieldOut,
        "x",
        "dvfOut.x"
    );
    eqns.evaluateDimensionedTypeField
    (
        dimensionedVectorFieldOut,
        "y",
        "dvfOut.y"
    );
    eqns.evaluateDimensionedTypeField
    (
        dimensionedVectorFieldOut,
        "z",
        "dvfOut.z"
    );

    eqns.evaluateGeometricTypeField(volVectorFieldOut, "x", "volVfOut.x");
    eqns.evaluateGeometricTypeField(volVectorFieldOut, "y", "volVfOut.y");
    eqns.evaluateGeometricTypeField(volVectorFieldOut, "z", "volVfOut.z");

    // tensors
    tensorOut.xx() = eqns.evaluateScalar("tOut.xx");
    tensorOut.xy() = eqns.evaluateScalar("tOut.xy");
    tensorOut.xz() = eqns.evaluateScalar("tOut.xz");
    tensorOut.yx() = eqns.evaluateScalar("tOut.yx");
    tensorOut.yy() = eqns.evaluateScalar("tOut.yy");
    tensorOut.yz() = eqns.evaluateScalar("tOut.yz");
    tensorOut.zx() = eqns.evaluateScalar("tOut.zx");
    tensorOut.zy() = eqns.evaluateScalar("tOut.zy");
    tensorOut.zz() = eqns.evaluateScalar("tOut.zz");

    dimensionedTensorOut.value().xx() = eqns.evaluateScalar("dtOut.xx");
    dimensionedTensorOut.value().xy() = eqns.evaluateScalar("dtOut.xy");
    dimensionedTensorOut.value().xz() = eqns.evaluateScalar("dtOut.xz");
    dimensionedTensorOut.value().yx() = eqns.evaluateScalar("dtOut.yx");
    dimensionedTensorOut.value().yy() = eqns.evaluateScalar("dtOut.yy");
    dimensionedTensorOut.value().yz() = eqns.evaluateScalar("dtOut.yz");
    dimensionedTensorOut.value().zx() = eqns.evaluateScalar("dtOut.zx");
    dimensionedTensorOut.value().zy() = eqns.evaluateScalar("dtOut.zy");
    dimensionedTensorOut.value().zz() = eqns.evaluateScalar("dtOut.zz");
    // There is currently no elegant way to check the dimensions of each
    // component of a "dimensionedType".  (dimensionedTypeFields and
    // GeometricTypeFields work automatically, though.)  This is how we do
    // it:
    dimensionedTensorOut.dimensions() = eqns.evaluateDimensions("dtOut.xx");
    dimensionedTensorOut.dimensions() = eqns.evaluateDimensions("dtOut.xy");
    dimensionedTensorOut.dimensions() = eqns.evaluateDimensions("dtOut.xz");
    dimensionedTensorOut.dimensions() = eqns.evaluateDimensions("dtOut.yx");
    dimensionedTensorOut.dimensions() = eqns.evaluateDimensions("dtOut.yy");
    dimensionedTensorOut.dimensions() = eqns.evaluateDimensions("dtOut.yz");
    dimensionedTensorOut.dimensions() = eqns.evaluateDimensions("dtOut.zx");
    dimensionedTensorOut.dimensions() = eqns.evaluateDimensions("dtOut.zy");
    dimensionedTensorOut.dimensions() = eqns.evaluateDimensions("dtOut.zz");

    eqns.evaluateDimensionedTypeField
    (
        dimensionedTensorFieldOut,
        "xx",
        "dtfOut.xx"
    );
    eqns.evaluateDimensionedTypeField
    (
        dimensionedTensorFieldOut,
        "xy",
        "dtfOut.xy"
    );
    eqns.evaluateDimensionedTypeField
    (
        dimensionedTensorFieldOut,
        "xz",
        "dtfOut.xz"
    );
    eqns.evaluateDimensionedTypeField
    (
        dimensionedTensorFieldOut,
        "yx",
        "dtfOut.yx"
    );
    eqns.evaluateDimensionedTypeField
    (
        dimensionedTensorFieldOut,
        "yy",
        "dtfOut.yy"
    );
    eqns.evaluateDimensionedTypeField
    (
        dimensionedTensorFieldOut,
        "yz",
        "dtfOut.yz"
    );
    eqns.evaluateDimensionedTypeField
    (
        dimensionedTensorFieldOut,
        "zx",
        "dtfOut.zx"
    );
    eqns.evaluateDimensionedTypeField
    (
        dimensionedTensorFieldOut,
        "zy",
        "dtfOut.zy"
    );
    eqns.evaluateDimensionedTypeField
    (
        dimensionedTensorFieldOut,
        "zz",
        "dtfOut.zz"
    );
    eqns.evaluateGeometricTypeField(volTensorFieldOut, "xx", "volTfOut.xx");
    eqns.evaluateGeometricTypeField(volTensorFieldOut, "xy", "volTfOut.xy");
    eqns.evaluateGeometricTypeField(volTensorFieldOut, "xz", "volTfOut.xz");
    eqns.evaluateGeometricTypeField(volTensorFieldOut, "yx", "volTfOut.yx");
    eqns.evaluateGeometricTypeField(volTensorFieldOut, "yy", "volTfOut.yy");
    eqns.evaluateGeometricTypeField(volTensorFieldOut, "yz", "volTfOut.yz");
    eqns.evaluateGeometricTypeField(volTensorFieldOut, "zx", "volTfOut.zx");
    eqns.evaluateGeometricTypeField(volTensorFieldOut, "zy", "volTfOut.zy");
    eqns.evaluateGeometricTypeField(volTensorFieldOut, "zz", "volTfOut.zz");

    // Update the equation output dictionary
    eqnOutputDict.set("sOut", scalarOut);
    eqnOutputDict.set("dsOut", dimensionedScalarOut);
    eqnOutputDict.set("vOut", vectorOut);
    eqnOutputDict.set("dvOut", dimensionedVectorOut);
    eqnOutputDict.set("tOut", tensorOut);
    eqnOutputDict.set("dtOut", dimensionedTensorOut);

