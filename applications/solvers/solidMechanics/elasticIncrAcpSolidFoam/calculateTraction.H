{
    surfaceVectorField n = mesh.Sf()/mesh.magSf();

    // traction = (n & fvc::interpolate(sigma));

    // surfaceTensorField sGradU =
    //     ((I - n*n) & fvc::interpolate(gradU));

    // traction =
    //     (2*mu + lambda)*snGradU
    //   - (mu + lambda)*(snGradU&(I - n*n))
    //   + mu*(sGradU & n)
    //   + lambda*tr(sGradU&(I - n*n))*n;

    // traction =
    //     (2*mu + lambda)*fvc::snGrad(U)
    //    - (mu + lambda)*(n & sGradU)
    //    + mu*(sGradU & n)
    //    + lambda*tr(sGradU)*n;

    // philipc
    // I am having trouble with back-calculation of interface tractions from solid interface
    // procedure (in multiMaterial.C), the tractions have quite large differences from each
    // side. Interpolating sigma is OK for now
    // traction = rheology.law().interfaceTraction(n, U, gradU, rheology.mu(), rheology.lambda());
#   include "calculateDEpsilonDSigma.H"
    traction = (n&fvc::interpolate(sigma+DSigma));

    // forAll(traction.boundaryField(), patchi)
    // {
    //     if (mesh.boundary()[patchi].type() == "cohesive")
    //     {
    //         forAll(traction.boundaryField()[patchi], facei)
    //         {
    //             Pout << "face " << facei << " with traction magnitude "
    //                 << mag(traction.boundaryField()[patchi][facei])/1e6 << " MPa and traction "
    //                 << traction.boundaryField()[patchi][facei]/1e6 << " MPa" << endl;
    //         }
    //     }
    // }
}
