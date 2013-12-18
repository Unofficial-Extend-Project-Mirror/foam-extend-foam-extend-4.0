Info << "Reading field Cs" << endl;
areaScalarField Cs
(
    IOobject
    (
        "Cs",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh
);

dimensioned<scalar> Cs0
(
    "Cs0",
    dimensionSet(1, -2, 0, 0, 0, 0, 0),
    1.0
);

const areaVectorField& R = aMesh.areaCentres();

Cs = Cs0*(1.0 + R.component(vector::X)/mag(R));


dimensioned<scalar> Ds
(
    "Ds",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    1.0
);


areaVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh,
    dimensioned<vector>("Us", dimVelocity, vector::zero)
);

dimensioned<scalar> Uinf("Uinf", dimVelocity, 1.0);

forAll (Us, faceI)
{
    Us[faceI].x() =
        Uinf.value()*(0.25*(3.0 + sqr(R[faceI].x()/mag(R[faceI]))) - 1.0);

    Us[faceI].y() =
        Uinf.value()*0.25*R[faceI].x()*R[faceI].y()/sqr(mag(R[faceI]));

    Us[faceI].z() =
        Uinf.value()*0.25*R[faceI].x()*R[faceI].z()/sqr(mag(R[faceI]));
}

Us -= aMesh.faceAreaNormals()*(aMesh.faceAreaNormals() & Us);


edgeScalarField phis
(
    IOobject
    (
        "phis",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearEdgeInterpolate(Us) & aMesh.Le()
);

