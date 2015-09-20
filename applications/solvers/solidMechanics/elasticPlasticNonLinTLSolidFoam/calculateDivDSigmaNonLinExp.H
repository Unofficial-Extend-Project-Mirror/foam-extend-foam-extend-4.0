if(divDSigmaNonLinExpMethod == "standard")
{
    divDSigmaNonLinExp = fvc::div
    (
        (
            mu*
            (
                (gradDU & gradU.T())
              + (gradU & gradDU.T())
              + (gradDU & gradDU.T())
            )
        )
      + (
            0.5*lambda*
            tr
            (
                (gradDU & gradU.T())
              + (gradU & gradDU.T())
              + (gradDU & gradDU.T())
            )*I
        )
      + (DSigma & gradU)
      + ((sigma + DSigma) & gradDU),
        "div(sigma)"
    );
}
else if(divDSigmaNonLinExpMethod == "surface")
{
     divDSigmaNonLinExp = fvc::div
     (
         mesh.magSf()*
         (
             (
                 muf*
                 (
                     n & fvc::interpolate
                     (
                         (gradU & gradDU.T())
                       + (gradDU & gradU.T())
                       + (gradDU & gradDU.T())
                     )
                 )
             )
           + (
                 0.5*lambdaf*
                 (
                     n*tr
                     (
                         fvc::interpolate
                         (
                             (gradU & gradDU.T())
                           + (gradDU & gradU.T())
                           + (gradDU & gradDU.T())
                         )
                     )
                 )
             )
           + (
                 n & fvc::interpolate
                 (
                     (DSigma & gradU)
                   + ((sigma + DSigma) & gradDU)
                 )
             )
         )
     );

     // divDSigmaNonLinExp =
     //   fvc::div
     //   (
     //     ( muf * mesh.Sf()
     //       & fvc::interpolate(
     //                  (gradDU & gradU.T())
     //                  + (gradU & gradDU.T())
     //                  + (gradDU & gradDU.T())
     //                  ) )
     //     + ( lambdaf * 0.5* tr(
     //                   fvc::interpolate(
     //                            (gradDU & gradU.T())
     //                            + (gradU & gradDU.T())
     //                            + (gradDU & gradDU.T())
     //                            )
     //                   ) * mesh.Sf() )
     //     + (mesh.Sf() & fvc::interpolate( (DSigma & gradU )
     //                      + ( (sigma + DSigma) & gradDU )
     //                      ) )
     //     );
}
else
{
       FatalErrorIn(args.executable())
           << "divDSigmaExp method " << divDSigmaExpMethod << " not found!"
           << abort(FatalError);
}
