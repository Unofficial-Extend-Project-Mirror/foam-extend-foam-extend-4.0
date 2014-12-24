//--------------------------------------------------//
//- rotate fields
//--------------------------------------------------//
{
  Info << "Rotating fields" << endl;

  volTensorField F = I + gradDU;

  volTensorField Finv = inv(F);

  volScalarField J = det(F);
  if(min(J.internalField()) < 0)
    {
      FatalErrorIn(args.executable())
        << "Negative Jacobian -> cell negative volumes!!!"
        << exit(FatalError);
    }

  rho = rho/J;
  n = mesh.Sf()/mesh.magSf();

  // rotate strain and stress fields
  //epsilon = symm(Finv & epsilon & Finv.T());
  epsilon = transform(Finv, epsilon);

  //sigma  = 1/J * symm(F.T() & sigma & F);
  sigma = (1/J)*transform(F.T(), sigma);

  // Note: for strict large strain approach, the constitutive stiffness
  // would be rotated here - see elasticOrthoNonLinULSolidFoam.
}
