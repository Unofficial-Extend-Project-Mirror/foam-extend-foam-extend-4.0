{
    Depsilon = symm(gradDU);

    DSigma = 2*mu*Depsilon + lambda*(I*tr(Depsilon));
}
