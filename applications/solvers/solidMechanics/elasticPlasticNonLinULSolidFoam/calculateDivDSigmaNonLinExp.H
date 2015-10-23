if(divDSigmaNonLinExpMethod == "standard")
{
    divDSigmaNonLinExp = fvc::div
    (
        ( mu * (gradDU & gradDU.T()) )
      + ( 0.5 * lambda * tr(gradDU & gradDU.T()) * I )
      + ( (sigma + DSigma) & gradDU ),
        "div(sigma)"
    );
}
else if(divDSigmaNonLinExpMethod == "surface")
{
    divDSigmaNonLinExp = fvc::div
    (
        mesh.magSf()*
        (
            ( muf * (n & fvc::interpolate( gradDU & gradDU.T() )) )
          + ( 0.5*lambdaf * (n * tr(fvc::interpolate( gradDU & gradDU.T() ))) )
          + (n & fvc::interpolate( (sigma + DSigma) & gradDU ))
        )
    );
}
else
{
    FatalError << "divDSigmaExp method " << divDSigmaExpMethod << " not found!" << endl;
}
