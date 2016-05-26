    Info<< "Reading thermophysical properties\n" << endl;

    volScalarField* pPtr = NULL;
    volScalarField* hPtr = NULL;
    const volScalarField* TPtr = NULL;

    volScalarField* psisPtr = NULL;
    basicThermo* thermoPtr = NULL;

    Switch compressible;

    {
        dictionary pimple = mesh.solutionDict().subDict("PIMPLE");

        compressible = Switch(pimple.lookup("compressible"));

        if (compressible)
        {
            thermoPtr =
            (
                basicPsiThermo::New(mesh)
            ).ptr();

            pPtr = &(thermoPtr->p());
            hPtr = &(thermoPtr->h());
            TPtr = &(thermoPtr->T());
            psisPtr = new volScalarField
            (
                "psi",
                thermoPtr->psi()/thermoPtr->Cp()*thermoPtr->Cv()
            );
        }
        else
        {
            thermoPtr =
            (
                basicRhoThermo::New(mesh)
            ).ptr();

            pPtr = &(thermoPtr->p());
            hPtr = &(thermoPtr->h());
            TPtr = &(thermoPtr->T());
            psisPtr = const_cast<volScalarField*>(&(thermoPtr->psi()));
        }
    }

    basicThermo& thermo = *thermoPtr;
    volScalarField& p = *pPtr;
    volScalarField& h = *hPtr;
    const volScalarField& T = *TPtr;
    volScalarField& psis = *psisPtr;
