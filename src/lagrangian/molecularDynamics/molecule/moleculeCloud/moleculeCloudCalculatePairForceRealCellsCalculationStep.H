idI = molI->id();

idJ = molJ->id();

rIJ = molI->position() - molJ->position();

rIJMagSq = magSqr(rIJ);

if(pairPotentials_.rCutSqr(idI, idJ, rIJMagSq))
{
    rIJMag = mag(rIJ);

    fIJ = (rIJ/rIJMag)*pairPotentials_.force(idI, idJ, rIJMag);

    // Acceleration increment for molI
    molI->A() += fIJ/(molI->mass());

    // Acceleration increment for molJ
    molJ->A() += -fIJ/(molJ->mass());

    scalar potentialEnergy
    (
        pairPotentials_.energy(idI, idJ, rIJMag)
    );

    molI->potentialEnergy() += 0.5*potentialEnergy;

    molJ->potentialEnergy() += 0.5*potentialEnergy;

    molI->rf() += rIJ * fIJ;

    molJ->rf() += rIJ * fIJ;
}
