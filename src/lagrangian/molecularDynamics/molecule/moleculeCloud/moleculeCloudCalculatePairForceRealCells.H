
vector rIJ;

scalar rIJMag;

scalar rIJMagSq;

vector fIJ;

label idI;

label idJ;

mol = this->begin();

molecule* molI = &mol();

molecule* molJ = &mol();

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

#               include "moleculeCloudCalculatePairForceRealCellsCalculationStep.H"
            }
        }

        forAll(cellOccupancy_[dIL],cellIOtherMols)
        {
            molJ = cellOccupancy_[dIL][cellIOtherMols];

            if (molJ > molI)
            {
#                include "moleculeCloudCalculatePairForceRealCellsCalculationStep.H"
            }
        }
    }
}
