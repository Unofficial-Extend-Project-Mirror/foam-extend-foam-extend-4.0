//--------------------------------------------------//
//- rotate fields
//--------------------------------------------------//
{
  Info << "Rotating Fields" << endl;

  //- dF is calculated at the end of the momentum loop
  F = I + DF;
  volTensorField Finv = inv(F);
  J = det(F);
  
  rho = rho/J;
  sigma = 1/J * symm(F & sigma & F.T());
  epsilon = symm(Finv.T() & epsilon & Finv);
  epsilonP = symm(Finv.T() & epsilonP & Finv);
  epsilonE = symm(Finv.T() & epsilonE & Finv);
}
