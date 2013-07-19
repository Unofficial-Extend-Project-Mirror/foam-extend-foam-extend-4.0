{
    // Solid side
#   include "readSolidControls.H"

    for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        coupledFvScalarMatrix TEqns(2);

        fvScalarMatrix* TFluidEqn = new fvScalarMatrix
        (
            rho*Cp*
            (
                fvm::div(phi, T)
              + fvm::SuSp(-fvc::div(phi), T)
            )
          - fvm::laplacian(kappaEff, T)
            ==
            radiation->Ru()
          - fvm::Sp(4.0*radiation->Rp()*pow3(T), T)
          + 3.0*radiation->Rp()*pow4(T)
        );

        //Called automatically in 1.6.1-ext
        //TFluidEqn->boundaryManipulate(T.boundaryField());
        TFluidEqn->relax();

        fvScalarMatrix* TSolidEqn = new fvScalarMatrix
        (
          - fvm::laplacian(ksolidf, Tsolid, "laplacian(k,T)")
          + fvm::SuSp(-solidThermo.S()/Tsolid, Tsolid)
        );

        //Called automatically in 1.6.1-ext
        //TSolidEqn->boundaryManipulate(Tsolid.boundaryField());
        TSolidEqn->relax();


        // Add fluid equation
        TEqns.set(0, TFluidEqn);

        // Add solid equation
        TEqns.set(1, TSolidEqn);

        TEqns.solve();
    }
}
