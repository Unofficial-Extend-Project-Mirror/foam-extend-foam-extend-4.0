{
    volSymmTensorField Depsilon = symm(gradDU)
      + 0.5*symm(gradDU & gradDU.T());

    DSigma = 2*mu*Depsilon + I*(lambda*tr(Depsilon));
}
