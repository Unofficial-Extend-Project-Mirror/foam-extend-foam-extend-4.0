//--------------------------------------------------//
//- rotate fields
//--------------------------------------------------//
{
  Info << "Rotating fields" << endl;

  F = I + gradDU;
  
  U += DU;
  
  epsilon += DEpsilon;

  sigma += DSigma;

  volTensorField Finv = inv(F);

  J = det(F);
  
  rho = rho/J;

  epsilon = symm(Finv & epsilon & Finv.T());

  sigma  = 1/J * symm(F.T() & sigma & F);
}
