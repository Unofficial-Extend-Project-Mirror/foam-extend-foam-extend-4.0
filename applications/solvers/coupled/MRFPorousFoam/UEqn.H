{
    // Update boundary velocity for consistency with the flux
    mrfZones.correctBoundaryVelocity(U);

    // Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    // Add MRF and porous sources implicitly.  HJ, 18/Nov/2017
    tmp<volTensorField> tTU;

    if (addMRF || addPorosity)
    {
        tTU = tmp<volTensorField>
        (
            new volTensorField
            (
                IOobject
                (
                    "TU",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedTensor("zero", dimless/dimTime, tensor::zero)
            )
        );
        volTensorField& TU = tTU();

        // Add implicit MRF source as a Hodge dual of the rotational velocity
        TU += *mrfZones.omega();
        
        // Add implicit resistance
        pZones.addResistance(UEqn, TU);

        trTU = inv(TU + tensor(I)*UEqn.A());
        trTU().rename("rAU");
    }
    else
    {
        trAU = 1.0/UEqn.A();
        trAU().rename("rAU");
    }

    // Under-relax momentum. Note this will destroy the H and A
    UEqn.relax();

    // Insert momentum equation
    UpEqn.insertEquation(0, UEqn);

    if (addMRF || addPorosity)
    {
        // Manually over-ride the 3x3 block to handle the off-diagonal
        // part of the Ap coefficient
        const tensorField& TUIn = tTU().internalField();

        CoeffField<vector4>::squareTypeField& DD = UpEqn.diag().asSquare();

        const scalarField& V = mesh.V().field();

        // Note: insertion should only happen in MRF and porous cell zones
        // HJ, 14/Mar/2016

        // Warning.  Possible problem with a zone which is both MRF and porous
        // The solution is to do the loop below everywhere
        // Currently only inserting zones for efficiency.  HJ, 18/Nov/2017
        register label cellI;

        forAll (mrfZones, mrfZoneI)
        {
            const labelList& curZoneCells = mrfZones[mrfZoneI].zone();

            // Loop over all cells in the zone
            forAll (curZoneCells, zcI)
            {
                cellI = curZoneCells[zcI];

                const scalar& cellV = V[cellI];

                const tensor& cellTU = TUIn[cellI];

                CoeffField<vector4>::squareType& cellDD = DD[cellI];

                cellDD(0, 0) += cellV*cellTU.xx();
                cellDD(0, 1) += cellV*cellTU.xy();
                cellDD(0, 2) += cellV*cellTU.xz();

                cellDD(1, 0) += cellV*cellTU.yx();
                cellDD(1, 1) += cellV*cellTU.yy();
                cellDD(2, 2) += cellV*cellTU.yz();

                cellDD(2, 0) += cellV*cellTU.zx();
                cellDD(2, 1) += cellV*cellTU.zy();
                cellDD(2, 2) += cellV*cellTU.zz();
            }
        }

        forAll (pZones, pZoneI)
        {
            const labelList& curZoneCells = pZones[pZoneI].zone();

            // Loop over all cells in the zone
            forAll (curZoneCells, zcI)
            {
                cellI = curZoneCells[zcI];

                const scalar& cellV = V[cellI];

                const tensor& cellTU = TUIn[cellI];

                CoeffField<vector4>::squareType& cellDD = DD[cellI];

                cellDD(0, 0) += cellV*cellTU.xx();
                cellDD(0, 1) += cellV*cellTU.xy();
                cellDD(0, 2) += cellV*cellTU.xz();

                cellDD(1, 0) += cellV*cellTU.yx();
                cellDD(1, 1) += cellV*cellTU.yy();
                cellDD(2, 2) += cellV*cellTU.yz();

                cellDD(2, 0) += cellV*cellTU.zx();
                cellDD(2, 1) += cellV*cellTU.zy();
                cellDD(2, 2) += cellV*cellTU.zz();
            }
        }
    }
}
