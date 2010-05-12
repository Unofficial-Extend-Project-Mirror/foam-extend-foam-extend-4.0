{ 
    vector rIJ;

    scalar rIJMag;

    scalar rIJMagSq;

    label idI;

    label idJ;

    mol = this->begin();

    molecule* molI = &mol();

    molecule* molJ = &mol();

    DynamicList<molecule*> molsToDelete;

    forAll(directInteractionList_, dIL)
    {
        forAll(cellOccupancy_[dIL],cellIMols)
        {
            molI = cellOccupancy_[dIL][cellIMols];

            forAll(directInteractionList_[dIL], interactingCells)
            {
                List< molecule* > cellJ =
                    cellOccupancy_[directInteractionList_[dIL][interactingCells]];

                forAll(cellJ, cellJMols)
                {
                    molJ = cellJ[cellJMols];

#                   include "moleculeCloudRemoveHighEnergyOverlapsRealCellsCalculationStep.H"
                }
            }

            forAll(cellOccupancy_[dIL],cellIOtherMols)
            {
                molJ = cellOccupancy_[dIL][cellIOtherMols];

                if (molJ > molI)
                {
#                   include "moleculeCloudRemoveHighEnergyOverlapsRealCellsCalculationStep.H"
                }
            }
        }
    }

    forAll (molsToDelete, mTD)
    {
        deleteParticle(*(molsToDelete[mTD]));
    }
}
