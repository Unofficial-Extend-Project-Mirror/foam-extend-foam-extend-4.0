if(divSigmaExpMethod == "standard")
{
    //- calculating the full gradient has good convergence and no high freq oscillations
    divSigmaExp = fvc::div((C && symm(gradU)) - (K & gradU), "div(sigma)");
}
else if(divSigmaExpMethod == "surface")
{
    //- this form seems to have the best convergence
    divSigmaExp = fvc::div
    (
        mesh.magSf()*
        (
            (n & (Cf && fvc::interpolate(symm(gradU))))
          - (n & (Kf & fvc::interpolate(gradU)))
        )
    );
}
else if(divSigmaExpMethod == "laplacian")
{
    //- can cause high freq oscillations and slow convergence
    divSigmaExp =
        fvc::div(C && symm(epsilon), "div(sigma)")
      - fvc::laplacian(K, U, "laplacian(K, U)");
}
else
{
    FatalErrorIn(args.executable())
        << "divSigmaExp method " << divSigmaExpMethod << " not found!" << endl;
}
