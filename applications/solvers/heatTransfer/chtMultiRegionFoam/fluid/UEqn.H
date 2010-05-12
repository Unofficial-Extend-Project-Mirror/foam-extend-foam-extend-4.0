    tmp<fvVectorMatrix> UEqn = solveMomentumEquation
    (
        momentumPredictor,
        Uf[i],
        rhof[i],
        phif[i],
        pdf[i],
        ghf[i],
        turb[i]
    );
