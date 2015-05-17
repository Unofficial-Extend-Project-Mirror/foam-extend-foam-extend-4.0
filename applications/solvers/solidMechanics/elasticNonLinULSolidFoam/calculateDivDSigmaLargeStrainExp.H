//----------------------------------------------------//
//- sigma explicit non linear explicit terms
//----------------------------------------------------//
if(divDSigmaLargeStrainExpMethod == "standard")
{
    divDSigmaLargeStrainExp = fvc::div
    (
        mu*(gradDU & gradDU.T())
        //+ 0.5*lambda*(gradDU && gradDU)*I //- equivalent to 0.5*lambda*(I*tr(gradDU & gradDU.T()))
      + 0.5*lambda*tr(gradDU & gradDU.T())*I
      + ((sigma + DSigma) & gradDU),
        "div(sigma)"
    );
}
else if(divDSigmaLargeStrainExpMethod == "surface")
{
    divDSigmaLargeStrainExp = fvc::div
    (
        muf * (mesh.Sf() & fvc::interpolate(gradDU & gradDU.T()))
      + 0.5*lambdaf * (mesh.Sf() & (fvc::interpolate(gradDU && gradDU)*I))
      + (mesh.Sf() & fvc::interpolate( (sigma + DSigma) & gradDU ))
    );
}
else
{
    FatalError
        << "divDSigmaLargeStrainMethod not found!"
        << exit(FatalError);
}
