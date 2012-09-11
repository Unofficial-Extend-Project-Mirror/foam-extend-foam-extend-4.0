DEpsilon = symm(gradDU) + 0.5*symm(gradDU & gradDU.T());

DSigma = 2*mu*DEpsilon + lambda*(I*tr(DEpsilon));
