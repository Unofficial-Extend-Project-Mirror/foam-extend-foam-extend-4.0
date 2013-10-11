//- Increment of Green finite strain tensor
DEpsilon = symm(gradDU) + 0.5*symm((gradDU & gradU.T()) + (gradU & gradDU.T()) + (gradDU & gradDU.T()));

//- Increment of second Piola-Kirchhoff stress tensor
DSigma = 2*mu*DEpsilon + lambda*(I*tr(DEpsilon));
