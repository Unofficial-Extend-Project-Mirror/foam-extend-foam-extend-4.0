if(Pstream::master())
{
    scalarField weights = 
        mag(phi.boundaryField()[interface.bPatchID()]);

    if(mag(sum(weights)) > SMALL)
    {
        weights /= sum(weights);
    }

    phi.boundaryField()[interface.bPatchID()] -=
        weights*sum(phi.boundaryField()[interface.bPatchID()]);
}
