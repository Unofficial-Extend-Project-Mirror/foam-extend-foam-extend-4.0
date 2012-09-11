DEpsilon = symm(gradDU) + 0.5*symm(gradDU & gradDU.T());
	
DSigma = 2*rheology.newMu()*(DEpsilon - DEpsilonP)
  + rheology.newLambda()*(I*tr(DEpsilon));
