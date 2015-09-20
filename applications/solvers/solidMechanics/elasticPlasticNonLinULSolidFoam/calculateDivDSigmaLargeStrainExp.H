//----------------------------------------------------//
//- sigma explicit large strain explicit terms
//----------------------------------------------------//
if(divDSigmaLargeStrainExpMethod == "standard")
{
    divDSigmaLargeStrainExp = fvc::div
    (
        mu*(gradDU & gradDU.T())
      + 0.5*lambda*(gradDU && gradDU)*I //- equivalent to 0.5*lambda*(I*tr(gradDU & gradDU.T()))
      + ((sigma + DSigma) & DF.T()),
        "div(sigma)"
    );
}
else if(divDSigmaLargeStrainExpMethod == "surface")
{
    divDSigmaLargeStrainExp = fvc::div
    (
        muf * (mesh.Sf() & fvc::interpolate(gradDU & gradDU.T()))
      + 0.5*lambdaf * (mesh.Sf() & (fvc::interpolate(gradDU && gradDU)*I))
      + (mesh.Sf() & fvc::interpolate( sigma & DF.T() ))
      + (mesh.Sf() & fvc::interpolate(DSigma & DF.T() ))
    );
}
else
{
    FatalError
        << "divDSigmaLargeStrainMethod not found!"
        << exit(FatalError);
}

//- relax
divDSigmaLargeStrainExp.relax();
