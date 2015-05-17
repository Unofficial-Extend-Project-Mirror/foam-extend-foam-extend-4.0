if(divDSigmaExpMethod == "standard")
{
    divDSigmaExp = fvc::div
    (
        mu*gradDU.T() + lambda*(I*tr(gradDU)) - (mu + lambda)*gradDU,
        "div(sigma)"
    );
}
else if(divDSigmaExpMethod == "surface")
{
    divDSigmaExp = fvc::div
    (
        muf*(mesh.Sf() & fvc::interpolate(gradDU.T()))
      + lambdaf*(mesh.Sf() & I*fvc::interpolate(tr(gradDU)))
      - (muf + lambdaf)*(mesh.Sf() & fvc::interpolate(gradDU))
    );
}
else if(divDSigmaExpMethod == "decompose")
{
    snGradDU = fvc::snGrad(DU);

    surfaceTensorField shearGradDU = ((I - n*n) & fvc::interpolate(gradDU));

    divDSigmaExp = fvc::div
    (
        mesh.magSf()*
        (
          - (muf + lambdaf)*(snGradDU & (I - n*n))
          + lambdaf*tr(shearGradDU & (I - n*n))*n
          + muf*(shearGradDU & n)
        )
    );
}
else if(divDSigmaExpMethod == "expLaplacian")
{
    divDSigmaExp =
      - fvc::laplacian(mu + lambda, DU, "laplacian(DDU,DU)")
      + fvc::div(mu*gradDU.T() + lambda*(I*tr(gradDU)), "div(sigma)");
}
else
{
    FatalErrorIn(args.executable())
        << "divDSigmaExp method " << divDSigmaExpMethod << " not found!"
        << abort(FatalError);
}
