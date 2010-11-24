
// Volume surfactant concentration

volScalarField* CPtr = NULL;

if(!interface.cleanInterface() && interface.surfactant().soluble())
{
    Info << "\nReading field C\n" << endl;
    CPtr = new volScalarField
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    if(runTime.timeIndex() == 0)
    {
        const dimensionedScalar& CA = 
            interface.surfactant().surfactBulkConc();

        dimensionedScalar CB("CB", CA.dimensions(), 0);

        (*CPtr) = fluidIndicator*(CA - CB) + CB;

        if
        (
            CPtr->boundaryField()[spacePatchID].type()
         == inletOutletFvPatchScalarField::typeName
        )
        {
            inletOutletFvPatchScalarField& spaceC =
                refCast<inletOutletFvPatchScalarField>
                (
                    CPtr->boundaryField()[spacePatchID]
                );

            spaceC.refValue() = CA.value();            
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Boundary condition for bulk surfactant concentration at "
                    << "space patch isn't and should be: "
                    << inletOutletFvPatchScalarField::typeName
                    << abort(FatalError);
        }
        
        CPtr->correctBoundaryConditions();

        CPtr->write();
    }
}
