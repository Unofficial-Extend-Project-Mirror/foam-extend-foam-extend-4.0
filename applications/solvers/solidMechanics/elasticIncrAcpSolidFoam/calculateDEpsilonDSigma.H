DEpsilon = symm(gradDU);

DSigma = 2*mu*DEpsilon + lambda*(I*tr(DEpsilon));
