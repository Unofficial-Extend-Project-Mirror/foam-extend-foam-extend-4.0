Solid Mechanics
Finite Volume Solvers

The included solid mechanics solvers employ the finite volume method
(not finite elements/elephants) to numerically approximate the
displacements and stresses in solid bodies undergoing deformation.

The included solvers feature the following capabilities:
    small strain
    small strain with large rotations
    large strain
    Mises-Levy J2 plasticity
    thermal-elasticity
    visco-elasticity
    gravity body forces
    fluid-structure interactions
    multi-material analyses
    contact stress analysis with friction
    small strain orthotropic elasticity
    large strain orthotropic elasticity
    cohesive zones
        predefined crack path
        arbitrary crack propagation
    custom boundary conditions
    Aitken's under-relaation for displacement field

A number of people have contributed to the development of the solvers,
mainly within Alojz Ivankovic's research group. The code has been
assembled and is maintained by Philip Cardiff (University College Dublin),
and significant contributions have been made by Aleksandar Karac, Zeljko
Tukovic, Hrvoje Jasak, Declan Carolan, Michael Leonard, Valentine
Kanyanta, David McAuliffe, Declan McNamara and Tian Tang.

Have fun.

Philip



The folowing references are relevant and citations are welcome:

Cardiff P, Karać A & Ivanković A, A Large Strain Finite Volume Method for
Orthotropic Bodies with General Material Orientations, Computer Methods
in Applied Mechanics & Engineering, 2013,
http://dx.doi.org/10.1016/j.cma.2013.09.008.

Cardiff P, Karać A & Ivanković A, Development of a finite volume contact
solver based on the penalty method. Computational Materials Science, 64
283-284, 2012, http://dx.doi.org/10.1016/j.commatsci.2012.03.011.

Cardiff P, Karać A, Tuković Z & Ivanković A, Development of a finite volume
based structural solver for large rotation of non-orthogonal meshes, 7th
OpenFOAM Workshop, Darmstadt, Germany, 2012.

Tuković Z, Ivanković A & Karać A, Finite volume stress analysis in multi-
material linear elastic body. International Journal for Numerical Methods
in Engineering, 2012. doi:10.1002/nme.

Carolan D, Tuković Z, Murphy N, Ivanković A, Arbitrary crack propagation
in multi-phase materials using the finite volume method, Computational
Materials Science, 2013, http://dx.doi.org/10.1016/j.commatsci.2012.11.049.

Tuković Z & Jasak H, Updated lagrangian finite volume solver for large
deformation dynamic response of elastic body. Transactions of FAMENA,
1(31):1–16, 2007.

Jasak H & Tuković Z, Dynamic mesh handling in OpenFOAM applied to fluid-
structure interaction simulations, 5th European Conference on Computational
Fluid Dynamics ECCOMAS CFD, Lisbon, Portugal, 2010.

Tuković Z & Jasak H, Finite volume method for fluid-strucutre-interaction
with large structural displacements, 2nd OpenFOAM Workshop, Zagreb, 2007.

Jasak H & Weller H, Finite volume methodology for contact problems of linear
elastic solids, 3rd International Conference of Croatian Society of Mechanics,
pages 253–260, Cavtat/Dubrovnik, Crotatia, 2000.

Jasak H & Weller H, Application of the finite volume method and unstructured
meshes to linear elasticity, International Journal for Numerical Methods in
Engineering, pages 267–287, 2000.

Maneeratana K, Development of the finite volume method for non-linear
structural applications, PhD thesis, Imperial College London, 2000.

Cardiff P, Development of the finite volume method for hip joint stress
analysis, PhD thesis, University College Dublin, 2012.

Tang T, Hededal O, Cardif P, Roenby J, A Finite Volume Method solver for
non-linear soil stress analysis using OpenFOAM, 8th OpenFOAM Workshop,
Jeju, 2013.
