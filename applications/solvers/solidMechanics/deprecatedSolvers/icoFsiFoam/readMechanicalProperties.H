    Info<< "Reading mechanical properties\n" << endl;

    IOdictionary mechanicalProperties
    (
        IOobject
        (
            "mechanicalProperties",
            runTime.constant(),
            stressMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rho(mechanicalProperties.lookup("rho"));
    dimensionedScalar rhoE(mechanicalProperties.lookup("E"));
    dimensionedScalar nuS(mechanicalProperties.lookup("nu"));

    Info<< "Normalising E : E/rho\n" << endl;
    dimensionedScalar E = rhoE/rho;

    Info<< "Calculating Lame's coefficients\n" << endl;

    dimensionedScalar mu = E/(2.0*(1.0 + nuS));
    dimensionedScalar lambda = nuS*E/((1.0 + nuS)*(1.0 - 2.0*nuS));

    Switch planeStress(mechanicalProperties.lookup("planeStress"));

    if (planeStress)
    {
        Info<< "Plane Stress\n" << endl;

        // Change lambda for plane stress
        lambda = nuS*E/((1.0 + nuS)*(1.0 - nuS));
    }
    else
    {
        Info<< "Plane Strain\n" << endl;
    }
    Info<< "mu = " << mu.value() << " Pa\n";
    Info<< "lambda = " << lambda.value() << " Pa\n";
