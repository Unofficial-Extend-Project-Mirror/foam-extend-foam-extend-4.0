    if (runTime.outputTime())
    {
        volScalarField ca = contact.contactArea();
        ca.write();
    }
