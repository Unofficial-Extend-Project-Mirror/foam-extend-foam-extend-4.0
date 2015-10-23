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

dimensionedScalar mug
(
    transportProperties.lookup("mug")
);

dimensionedScalar mul
(
    transportProperties.lookup("mul")
);

dimensionedScalar sigma
(
    transportProperties.lookup("sigma")
);

dimensionedScalar rhol
(
    transportProperties.lookup("rhol")
);

dimensionedScalar rhog
(
    transportProperties.lookup("rhog")
);

dimensionedScalar h0
(
    transportProperties.lookup("h0")
);
