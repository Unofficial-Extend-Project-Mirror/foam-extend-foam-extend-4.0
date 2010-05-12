idI = molI->id();

idJ = molJ->id();

rIJ = molI->position() - molJ->position();

rIJMagSq = magSqr(rIJ);

if (pairPotentials_.rCutSqr(idI, idJ, rIJMagSq))
{
    rIJMag = mag(rIJ);

    bool remove = false;

    // Guard against pairPotentials_.energy being evaluated
    // if rIJMag < SMALL. A floating point exception will
    // happen otherwise.

    if (rIJMag < SMALL)
    {
        WarningIn("moleculeCloud::removeHighEnergyOverlaps()")
            << "Real molecule pair "
            << " idI = " << idI
            << " at position " << molI->position()
            << " idJ = " << idJ
            << " at position " << molJ->position()
            << " are closer than " << SMALL
            << ": mag separation = " << rIJMag
            << ". These may have been placed on top of each"
            << " other by a rounding error in molConfig in parallel"
            << " or a block filled with molecules twice."
            << " Removing one of the molecules."
            << endl;

        remove = true;
    }
    
    // Guard against pairPotentials_.energy being evaluated
    // if rIJMag < rMin. A tabulation lookup error will occur otherwise.

    if (rIJMag < pairPotentials_.rMin(idI, idJ))
    {
        remove = true;
    }

    if (!remove)
    {
        if
        (
            pairPotentials_.energy(idI, idJ, rIJMag) > potentialEnergyLimit_
        )
        {
           
            remove = true;
        }
    }

    if (remove)
    {
        if
        (
            idI == idJ
         || findIndex(removalOrder_, idJ) < findIndex(removalOrder_, idI)
        )
        {
            if (findIndex(molsToDelete, molJ) == -1)
            {
                molsToDelete.append(molJ);
            }
        }
        else
        {
            if (findIndex(molsToDelete, molI) == -1)
            {
                molsToDelete.append(molI);
            }
        }
    }
}
