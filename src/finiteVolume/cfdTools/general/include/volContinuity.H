if (mesh.moving())
{
    // The ddt term constructed by hand because it must be Euler
    const scalar dt = runTime.deltaT().value();

    const scalarField& V = mesh.V().field();
    const scalarField& V0 = mesh.V0().field();
    scalar sumV = gSum(V);
    scalar sumV0 = gSum(V0);

    // Rewrite of mesh motion check based on round-off error problems
    // Note precision around V and avoiding division with dt
    // HJ, 4/Oct/2011
    scalarField volChange = (V - V0)/V;

    scalarField divFlux =
      - fvc::div(mesh.phi())().internalField()*dt;

    scalarField conserve = volChange + divFlux;

    scalar sumLocalContErr = dt*gSum(mag(conserve*mesh.V()))/sumV;
    scalar globalContErr = dt*gSum(conserve*mesh.V())/sumV;

    label nMotionErrors = 0;
    scalar maxMotionError = 0;

    forAll (conserve, cellI)
    {
        maxMotionError = Foam::max(maxMotionError, mag(conserve[cellI]));

        if (mag(conserve[cellI]) > 1e-8)
        {
            Info<< "Motion conservation error in cell " << cellI << ": "
                << conserve[cellI]
                << " V: " << V[cellI] << " V0: " << V0[cellI]
                << " volume change: " << volChange[cellI]
                << " div flux: " << divFlux[cellI] << endl;

            nMotionErrors++;
        }
    }

    if (nMotionErrors > 0)
    {
        Info<< "Number of motion errors: " << nMotionErrors
            << " out of nCells = " << mesh.nCells() << endl;
    }

    Info<< "volume continuity errors : "
        << "volume = " << sumV
        << ", old volume = " << sumV0
        << ", max error = " << maxMotionError
        << ", sum local = " << sumLocalContErr
        << ", global = " << globalContErr << endl;
}
