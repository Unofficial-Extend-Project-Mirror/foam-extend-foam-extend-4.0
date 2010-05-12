if (runTime.outputTime())
{
    // Local free surface curvature and flux
    if(Pstream::master())
    {
        OFstream file
            (
                runTime.timePath()/"fsData.dat"
            );
        
        file.precision(12);

        scalarField flowRate =
            phi.boundaryField()[interface.aPatchID()]
          + fvc::meshPhi(rho,U)().boundaryField()[interface.aPatchID()];

        const vectorField& C = 
            mesh.Cf().boundaryField()[interface.aPatchID()];

        const scalarField& K = 
            interface.aMesh().faceCurvatures().internalField();

        file<< "X" << tab 
            << "Y" << tab 
            << "Z" << tab 
            << "K" << tab 
            << "flux";
            
        if(!interface.cleanInterface())
        {
            file<< tab << "Cs" << tab
                << "surfTension" << endl;
        }
        else
        {
            file << endl;
        }


        forAll(flowRate, faceI)
        {
            file<< C[faceI].x() << tab 
                << C[faceI].y() << tab
                << C[faceI].z() << tab
                << K[faceI] << tab
                << flowRate[faceI];

            if(!interface.cleanInterface())
            {
                file<< tab << 
                    interface.surfactantConcentration().internalField()[faceI] 
                    << tab
                    << interface.surfaceTension().internalField()[faceI] 
                    << endl;
            }
            else
            {
                file << endl;
            }
        }
    }
}
