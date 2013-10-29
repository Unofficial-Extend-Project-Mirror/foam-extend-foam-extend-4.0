DEpsilon = symm(gradDU);

DSigma = 2*mu*(DEpsilon - DEpsilonP) + lambda*(I*tr(DEpsilon));
