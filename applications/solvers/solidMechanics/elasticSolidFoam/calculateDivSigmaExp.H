if(divSigmaExpMethod == "standard")
{
    divSigmaExp = fvc::div
    (
        mu*gradU.T() + lambda*(I*tr(gradU)) - (mu + lambda)*gradU,
        "div(sigma)"
    );
}
else if(divSigmaExpMethod == "surface")
{
    divSigmaExp = fvc::div
    (
        muf*(mesh.Sf() & fvc::interpolate(gradU.T()))
      + lambdaf*(mesh.Sf() & I*fvc::interpolate(tr(gradU)))
      - (muf + lambdaf)*(mesh.Sf() & fvc::interpolate(gradU))
    );
}
else if(divSigmaExpMethod == "decompose")
{
    snGradU = fvc::snGrad(U);

    surfaceTensorField shearGradU = ((I - n*n) & fvc::interpolate(gradU));

    divSigmaExp = fvc::div
    (
        mesh.magSf()*
        (
          - (muf + lambdaf)*(snGradU & (I - n*n))
          + lambdaf*tr(shearGradU & (I - n*n))*n
          + muf*(shearGradU & n)
        )
    );
}
else if(divSigmaExpMethod == "expLaplacian")
{
    divSigmaExp =
      - fvc::laplacian(mu + lambda, U, "laplacian(DU,U)")
      + fvc::div(mu*gradU.T() + lambda*(I*tr(gradU)), "div(sigma)");
}
else
{
    FatalErrorIn(args.executable())
        << "divSigmaExp method " << divSigmaExpMethod << " not found!" << endl;
}
