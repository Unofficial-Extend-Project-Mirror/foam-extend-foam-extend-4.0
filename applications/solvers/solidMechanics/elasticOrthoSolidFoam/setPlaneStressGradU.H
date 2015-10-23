//- set gradU.zz() for plane stress
if(rheology.planeStress())
{
    forAll(gradU.internalField(), celli)
    {
        gradU.internalField()[celli].zz() =
            (
              - C.internalField()[celli].xxzz()*epsilon.internalField()[celli].xx()
              - C.internalField()[celli].yyzz()*epsilon.internalField()[celli].yy()
            )
            /
            C.internalField()[celli].zzzz();
    }
    gradU.correctBoundaryConditions();
}
