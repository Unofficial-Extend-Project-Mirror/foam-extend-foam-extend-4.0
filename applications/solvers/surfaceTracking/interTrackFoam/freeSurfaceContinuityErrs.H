if(Pstream::master())
{
    Info<< "Free surface flux: sum local = "
        << sum(mag(phi.boundaryField()[interface.aPatchID()])) 
        << ", global = "
        << sum(phi.boundaryField()[interface.aPatchID()]) << endl;

    Info<< "Free surface continuity error : sum local = " <<
        sum
        (
            mag
            (
                phiNet.boundaryField()[interface.aPatchID()]/
                interface.rhoFluidA().value()
            )
        ) 
        << ", global = " <<
        sum
        (
            phiNet.boundaryField()[interface.aPatchID()]/
            interface.rhoFluidA().value()
        ) 
        << endl << endl;
}
