//- Green finite strain tensor
epsilon = symm(gradU) + 0.5*symm(gradU & gradU.T());

//- second Piola-Kirchhoff stress tensor
sigma = 2*mu*epsilon + lambda*(I*tr(epsilon));
