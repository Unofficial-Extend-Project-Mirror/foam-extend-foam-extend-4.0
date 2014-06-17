{
    // Set reference
    for (label regI = 0; regI < mesh.regions().nRegions(); regI++)
    {
        if (needRef[regI])
        {
            UEqn.source()[refCells[regI]] +=
                UEqn.diag()[refCells[regI]]*vector::zero;

            UEqn.diag()[refCells[regI]] +=
                UEqn.diag()[refCells[regI]];
        }
    }

    // Set value for dead cells (small particles or regions)
//     label size = mesh.nCells() - label(sum(gamma.internalField()));

    if (deadCells.size())
    {
        labelList curRefCell(deadCells.size(), -1);
        label index = 0;
        forAll(deadCells, cellI)
        {
            curRefCell[index++] = deadCells[cellI];

//             if (gamma.internalField()[cellI] < 1 - SMALL)
//             {
//                 curRefCell[index++] = cellI;
//             }
        }

        vectorField curRefValue(curRefCell.size(), vector::zero);

        UEqn.setValues(curRefCell, curRefValue);
    }
}
