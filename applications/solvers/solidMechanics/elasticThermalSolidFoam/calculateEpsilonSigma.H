epsilon = symm(gradU);

sigma = 2*mu*epsilon + lambda*(I*tr(epsilon)) - I*threeKalpha*(T-T0);
