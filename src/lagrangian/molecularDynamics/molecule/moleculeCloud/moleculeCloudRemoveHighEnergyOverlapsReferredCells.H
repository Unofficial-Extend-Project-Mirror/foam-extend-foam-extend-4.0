{
    vector rKL;

    scalar rKLMag;

    scalar rKLMagSq;

    label idK;

    label idL;

    molecule* molK = &mol();

    DynamicList<molecule*> molsToDelete;

    forAll(referredInteractionList_, rIL)
    {
        referredCell& refCell = referredInteractionList_[rIL];

        forAll(refCell, refMols)
        {
            referredMolecule* molL = &(refCell[refMols]);

            List <label> realCells = refCell.realCellsForInteraction();

            forAll(realCells, rC)
            {
                label cellK = realCells[rC];

                List<molecule*> cellKMols = cellOccupancy_[cellK];

                forAll(cellKMols, cKM)
                {
                    molK = cellKMols[cKM];

                    idK = molK->id();

                    idL = molL->id();

                    rKL = molK->position() - molL->position();

                    rKLMagSq = magSqr(rKL);

                    if (pairPotentials_.rCutSqr(idK, idL, rKLMagSq))
                    {
                        rKLMag = mag(rKL);

                        bool remove = false;

                        // Guard against pairPotentials_.energy being evaluated
                        // if rKLMag < SMALL. A floating point exception will
                        // happen otherwise.

                        if (rKLMag < SMALL)
                        {
                            WarningIn
                            (
                                "moleculeCloud::removeHighEnergyOverlaps()"
                            )
                                << "Real-referred molecule pair "
                                << " idK = " << idK << " (real)"
                                << " at position " << molK->position()
                                << " idL = " << idL << " (referred)"
                                << " at position " << molL->position()
                                << " are closer than " << SMALL
                                << ": mag separation = " << rKLMag
                                << ". These may have been placed on top of each"
                                << " other by a rounding error in molConfig in"
                                << " parallel, a block filled with molecules"
                                << " twice, or a lattice ending very close"
                                << " to the edge of the mesh."
                                << " Removing one of the molecules."
                                << endl;

                                remove = true;
                        }

                        // Guard against pairPotentials_.energy being evaluated
                        // if rIJMag < rMin. A tubulation lookup error will occur otherwise.

                        if (rKLMag < pairPotentials_.rMin(idK, idL))
                        {
                            remove = true;
                        }

                        if (!remove)
                        {
                            if
                            (
                                pairPotentials_.energy(idK, idL, rKLMag)
                              > potentialEnergyLimit_
                            )
                            {
                                remove = true;
                            }
                        }

                        if (remove)
                        {
                            if
                            (
                                findIndex(removalOrder_, idK)
                              < findIndex(removalOrder_, idL)
                            )
                            {
                                if (findIndex(molsToDelete, molK) == -1)
                                {
                                    molsToDelete.append(molK);
                                }
                            }

                            else if
                            (
                                findIndex(removalOrder_, idK)
                             == findIndex(removalOrder_, idL)
                            )
                            {
                                if
                                (
                                    Pstream::myProcNo() == refCell.sourceProc()
                                 && cellK <= refCell.sourceCell()
                                )
                                {
                                    if (findIndex(molsToDelete, molK) == -1)
                                    {
                                        molsToDelete.append(molK);
                                    }
                                }

                                else if
                                (
                                    Pstream::myProcNo() < refCell.sourceProc()
                                )
                                {
                                    if (findIndex(molsToDelete, molK) == -1)
                                    {
                                        molsToDelete.append(molK);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    forAll (molsToDelete, mTD)
    {
        deleteParticle(*(molsToDelete[mTD]));
    }
}
