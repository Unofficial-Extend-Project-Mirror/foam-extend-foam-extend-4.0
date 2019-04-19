{
    while (simpleSolid.correctNonOrthogonal())
    {
        coupledFvScalarMatrix TEqns(2);

        // Fluid side
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

        TFluidEqn->relax();

        // Solid side
        fvScalarMatrix* TSolidEqn = new fvScalarMatrix
        (
          - fvm::laplacian(kSolidf, Tsolid, "laplacian(k,T)")
          + fvm::SuSp(-solidThermo.S()/Tsolid, Tsolid)
        );

        TSolidEqn->relax();

        // Add fluid equation
        TEqns.set(0, TFluidEqn);

        // Add solid equation
        TEqns.set(1, TSolidEqn);

        TEqns.solve();
    }
}
