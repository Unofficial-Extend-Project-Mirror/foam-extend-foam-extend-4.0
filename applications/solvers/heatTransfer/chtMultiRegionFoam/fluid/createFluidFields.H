    // Initialise fluid field pointer lists
    PtrList<basicThermo> thermof(fluidRegions.size());
    PtrList<volScalarField> rhof(fluidRegions.size());
    PtrList<volScalarField> Kf(fluidRegions.size());
    PtrList<volVectorField> Uf(fluidRegions.size());
    PtrList<surfaceScalarField> phif(fluidRegions.size());
    PtrList<compressible::RASModel> turb(fluidRegions.size());
    PtrList<volScalarField> DpDtf(fluidRegions.size());
    PtrList<volScalarField> ghf(fluidRegions.size());
    PtrList<volScalarField> pdf(fluidRegions.size());

    List<scalar> initialMassf(fluidRegions.size());

    dimensionedScalar pRef
    (
        "pRef",
        dimensionSet(1, -1, -2, 0, 0),
        rp.lookup("pRef")
    );

    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to pdf\n" << endl;
        pdf.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "pd",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to thermof\n" << endl;

        thermof.set
        (
            i,
            basicThermo::New(fluidRegions[i]).ptr()
        );

        Info<< "    Adding to rhof\n" << endl;
        rhof.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                thermof[i].rho()
            )
        );

        Info<< "    Adding to Kf\n" << endl;
        Kf.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "K",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermof[i].Cp()*thermof[i].alpha()
            )
        );

        Info<< "    Adding to Uf\n" << endl;
        Uf.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to phif\n" << endl;
        phif.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(rhof[i]*Uf[i])
                    & fluidRegions[i].Sf()
            )
        );

        Info<< "    Adding to turb\n" << endl;
        turb.set
        (
            i,
            autoPtr<compressible::RASModel>
            (
                compressible::RASModel::New
                (
                    rhof[i],
                    Uf[i],
                    phif[i],
                    thermof[i]
                )
            ).ptr()
        );

        Info<< "    Adding to DpDtf\n" << endl;
        DpDtf.set
        (
            i,
            new volScalarField
            (
                fvc::DDt
                (
                    surfaceScalarField
                    (
                        "phiU",
                        phif[i]/fvc::interpolate(rhof[i])
                    ),
                    thermof[i].p()
                )
            )
        );

        const dictionary& environmentalProperties =
            fluidRegions[i].lookupObject<IOdictionary>
                ("environmentalProperties");
        dimensionedVector g(environmentalProperties.lookup("g"));

        Info<< "    Adding to ghf\n" << endl;
        ghf.set
        (
            i,
            new volScalarField
            (
                "gh",
                 g & fluidRegions[i].C()
            )
        );

        Info<< "    Updating p from pd\n" << endl;
        thermof[i].p() == pdf[i] + rhof[i]*ghf[i] + pRef;
        thermof[i].correct();

        initialMassf[i] = fvc::domainIntegrate(rhof[i]).value();
    }


