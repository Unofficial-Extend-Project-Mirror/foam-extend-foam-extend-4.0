{
  // force residual is the net force on the model
  // this should got to zero in a converged steady state model
  // should be altered for parallel runs
  vector netForce = vector::zero;
  forAll(mesh.boundary(), patchi)
    {
      netForce +=
    sum(
        mesh.Sf().boundaryField()[patchi]
        &
        (
         2*mu.boundaryField()[patchi]*symm(gradU.boundaryField()[patchi])
         + lambda*tr(gradU.boundaryField()[patchi])*I
         )
        );
    }
  forceResidual = mag(netForce);
}
