if(interface.twoFluids() && pNeedRef)
{
    phi.boundaryField()[interface.bPatchID()] =
        (
            U.boundaryField()[interface.bPatchID()] 
          & mesh.Sf().boundaryField()[interface.bPatchID()]
        );
        
    scalarField weights = 
        mag(phi.boundaryField()[interface.bPatchID()]);

    if(mag(gSum(weights)) > VSMALL)
    {
        weights /= gSum(weights);
    }

    phi.boundaryField()[interface.bPatchID()] -=
        weights*gSum(phi.boundaryField()[interface.bPatchID()]);

//     Info << "Free surface shadow flow rate: "
//         << gSum(phi.boundaryField()[interface.bPatchID()]) << endl;

    phi.boundaryField()[interface.bPatchID()] +=
        p.boundaryField()[interface.bPatchID()].snGrad()
       *mesh.magSf().boundaryField()[interface.bPatchID()]
       /AU.boundaryField()[interface.bPatchID()];
}

