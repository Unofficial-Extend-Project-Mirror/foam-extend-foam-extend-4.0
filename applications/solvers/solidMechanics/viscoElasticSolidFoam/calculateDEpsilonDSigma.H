{
    DEpsilon = symm(gradDU);

    DSigma = 2.0*mu*DEpsilon  + lambda*(I*tr(DEpsilon));
}
