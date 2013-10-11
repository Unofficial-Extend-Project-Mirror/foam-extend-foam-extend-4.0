Info << "\nTime = " << runTime.timeName() 
    << ", iteration: " << outerCorr << endl;

if (outerCorr < 3 || couplingScheme == "FixedRelaxation")
{
    Info << "Current fsi under-relaxation factor: " 
        << fsiRelaxationFactor << endl;

    fluidPatchPointsDisplOld = fluidPatchPointsDispl;

    fluidPatchPointsDispl += fsiRelaxationFactor*fsiResidual;
}
else
{
    if (couplingScheme == "Aitken")
    {
        fsiRelaxationFactor = 
           -fsiRelaxationFactor
           *(
                gSum
                (
                    fsiResidualOld
                   &(fsiResidual - fsiResidualOld)
                )
               /(
                    gSum
                    (
                        (fsiResidual - fsiResidualOld)
                       &(fsiResidual - fsiResidualOld)
                    )
                )
           );

        fsiRelaxationFactor = mag(fsiRelaxationFactor);

        Info << "Current fsi under-relaxation factor (Aitken): " 
            << fsiRelaxationFactor << endl;        

        fluidPatchPointsDisplOld = fluidPatchPointsDispl;

        fluidPatchPointsDispl += 
            fsiRelaxationFactor*fsiResidual;
    }
//     else if (couplingScheme == "IQN-ILS")
//     {

//     }
}
