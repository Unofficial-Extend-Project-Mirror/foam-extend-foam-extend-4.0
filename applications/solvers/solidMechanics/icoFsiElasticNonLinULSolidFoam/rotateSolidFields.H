//--------------------------------------------------//
//- rotate fields
//--------------------------------------------------//
{
  Info << "Rotating fields" << endl;

  F = I + DF;
  
  //U += DU;
  Usolid += DU;
  
  //epsilon += DEpsilon;

  sigma += DSigma;

  //volTensorField Finv = inv(F);

  J = det(F);
  
  rho = rho/J;
 
  //epsilon = symm(Finv.T() & epsilon & Finv);

  sigma  = 1/J * symm(F & sigma & F.T());
}
