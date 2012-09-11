DEpsilon = symm(gradDU);

DSigma = 2*rheology.newMu()*(DEpsilon - DEpsilonP)
  + rheology.newLambda()*(I*tr(DEpsilon));
