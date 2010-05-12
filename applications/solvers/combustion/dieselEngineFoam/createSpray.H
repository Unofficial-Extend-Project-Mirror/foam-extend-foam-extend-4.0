Info << "Constructing Spray" << endl;

pointMesh pMesh(mesh);
volPointInterpolation vpi(mesh, pMesh);

PtrList<specieProperties> gasProperties(Y.size());
forAll(gasProperties, i)
{
    gasProperties.set
    (
        i,
        new specieProperties
        (
            dynamic_cast<const reactingMixture&>(thermo()).speciesData()[i]
        )
    );
}

spray dieselSpray
(
    vpi,
    U,
    rho,
    p,
    T,
    composition,
    gasProperties,
    thermo,
    environmentalProperties
);

scalar gasMass0 = fvc::domainIntegrate(rho).value();

if (dieselSpray.twoD())
{
    gasMass0 *= 2.0*mathematicalConstant::pi/dieselSpray.angleOfWedge();
}

gasMass0 -=
    dieselSpray.injectedMass(runTime.value()) - dieselSpray.liquidMass();
