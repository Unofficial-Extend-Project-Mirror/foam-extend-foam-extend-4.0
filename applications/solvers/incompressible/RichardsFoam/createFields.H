Info<< "Reading transportProperties\n" << endl;

// reading of the physical constant associated with the considered problem.
// Localisation of the data within the considered case:
// constant/transportProperties
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar K
(
    transportProperties.lookup("K")
);

dimensionedScalar alpha
(
    transportProperties.lookup("alpha")
);

dimensionedScalar thetas
(
    transportProperties.lookup("thetas")
);

dimensionedScalar thetar
(
    transportProperties.lookup("thetar")
);

dimensionedScalar n
(
    transportProperties.lookup("n")
);

dimensionedScalar C
(
    transportProperties.lookup("C")
);

dimensionedScalar S
(
    transportProperties.lookup("S")
);

// declaration of the variable and results fields

// Water velocity field [m/s]
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Water saturation field [-]
Info<< "Reading field theta\n" << endl;
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Water pressure field [m] - field of resolution.
Info<< "Reading field psi\n" << endl;
volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Field of residuals for the Picard loop [m].
Info<< "Reading field err\n" << endl;
volScalarField err
(
    IOobject
    (
        "err",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Dimensionned unit scalar field [m].
Info<< "Reading field vuz\n" << endl;
volVectorField vuz
(
    IOobject
    (
        "vuz",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Dimensionless unit vertical upward vector field.
Info<< "Reading field usf\n" << endl;
volScalarField usf
(
    IOobject
    (
        "usf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#   include "createPhi.H"

// Initialisation of the scalar containing the number of mesh cells. Note the
// use of gSum instead of sum.
double nbMesh;
nbMesh = gSum(usf);

// Initialisation of the scalar containing the residual for the exit test of the
// Picard loop.
double crit;
crit=0.;

// Initialisation of the token which counts the number of Picard iteraion for
// the adaptive time step procedure.
int currentPicard;
currentPicard = nIterPicard-3;

// Initialisation of the token which counts the number of Stabilisation cycle
// for the stabilisation of the adaptive time step procedure.
int sc;
sc = 0;

// Initialisation of the field of altitudes.
volVectorField positionVector = mesh.C();
volScalarField z = positionVector.component(vector::Z);

// Initialisation of the intermediate fields for the Picard loop.
volScalarField psi_tmp = psi;
volScalarField psim1 = psi;

// Initialisation of the varying transport properties for the Picard loop.
volScalarField thtil =
    0.5*
    (
        (1 + sign(psi)) + (1 - sign(psi))*
        pow((1 + pow(mag(alpha*psi),n)), - (1 - (1/n)))
    );

volScalarField thtil_tmp =
    0.5*
    (
        (1 + sign(psi_tmp)) + (1-sign(psi_tmp))*
        pow((1 + pow(mag(alpha*psi_tmp),n)), - (1 - (1/n)))
    );

volScalarField Krel =
    0.5*
    (
        (1 + sign(psi))*K + (1 - sign(psi))*K*pow(thtil,0.5)*
        pow((1 - pow((1 - pow(thtil,(n/(n - 1)))),(1 - (1/n)))),2)
    );

volScalarField Crel =
    S + 0.5*
    (
        (1 - sign(psi))*((thetas - thetar)*(thtil - thtil_tmp)*
        (1./((usf*pos(psi - psi_tmp)*pos(psi_tmp - psi)) + psi - psi_tmp)))
    );

// Initialisation of the gravity term.
volVectorField gradk = fvc::grad(Krel);
volScalarField gradkz = gradk.component(vector::Z);

// Initialisation of the velocity field.
U = - Krel*((fvc::grad(psi)) + vuz);
