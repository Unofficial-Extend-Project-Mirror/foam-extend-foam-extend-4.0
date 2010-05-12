vector rKL;

scalar rKLMag;

scalar rKLMagSq;

vector fKL;

label idK;

label idL;

molecule* molK = &mol();

forAll(referredInteractionList_, rIL)
{
    const List<label>& realCells =
        referredInteractionList_[rIL].realCellsForInteraction();

    forAll(referredInteractionList_[rIL], refMols)
    {
        referredMolecule* molL = &(referredInteractionList_[rIL][refMols]);

        forAll(realCells, rC)
        {
            List<molecule*> cellK = cellOccupancy_[realCells[rC]];

            forAll(cellK, cellKMols)
            {
                molK = cellK[cellKMols];

                idK = molK->id();

                idL = molL->id();

                rKL = molK->position() - molL->position();

                rKLMagSq = magSqr(rKL);

                if (pairPotentials_.rCutSqr(idK, idL, rKLMagSq))
                {
                    rKLMag = mag(rKL);

                    fKL = (rKL/rKLMag)*pairPotentials_.force(idK, idL, rKLMag);

                    // Acceleration increment for molK
                    molK->A() += fKL/(molK->mass());

                    // Adding a contribution of 1/2 of the potential energy
                    // from this interaction
                    molK->potentialEnergy() +=
                        0.5*pairPotentials_.energy(idK, idL, rKLMag);

                    molK->rf() += rKL * fKL;
                }
            }
        }
    }
}
