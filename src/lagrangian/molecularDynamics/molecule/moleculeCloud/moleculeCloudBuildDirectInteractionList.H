Info << nl << "Building list of direct interaction neighbours" << endl;

forAll (mesh_.points(), p)
{
    forAll(mesh_.faces(), f)
    {
        if(testPointFaceDistance(p, f))
        {
            const labelList& pCells(mesh_.pointCells()[p]);

            const label cellO(mesh_.faceOwner()[f]);

            forAll(pCells, pC)
            {
                const label cellI(pCells[pC]);

                // cells are not added to their own DIL

                if (cellO > cellI)
                {
                    if (findIndex(directInteractionList[cellI], cellO) == -1)
                    {
                        directInteractionList[cellI].append(cellO);
                    }
                }

                if (cellI > cellO)
                {
                    if (findIndex(directInteractionList[cellO], cellI) == -1)
                    {
                        directInteractionList[cellO].append(cellI);
                    }
                }

                if (mesh_.isInternalFace(f))
                {
                    // boundary faces will not have neighbour information
                    const label cellN(mesh_.faceNeighbour()[f]);

                    if (cellN > cellI)
                    {
                        if
                        (
                             findIndex(directInteractionList[cellI], cellN)
                          == -1
                        )
                        {
                            directInteractionList[cellI].append(cellN);
                        }
                    }

                    if (cellI > cellN)
                    {
                        if
                        (
                            findIndex(directInteractionList[cellN], cellI)
                         == -1
                        )
                        {
                            directInteractionList[cellN].append(cellI);
                        }
                    }
                }
            }
        }
    }
}

label edgeJIndex;

forAll (mesh_.edges(), edgeIIndex)
{
    const edge& eI(mesh_.edges()[edgeIIndex]);

    for
    (
        edgeJIndex = edgeIIndex + 1;
        edgeJIndex != mesh_.edges().size();
        ++edgeJIndex
    )
    {
        const edge& eJ(mesh_.edges()[edgeJIndex]);

        if (testEdgeEdgeDistance(eI, eJ))
        {
            const labelList& eICells(mesh_.edgeCells()[edgeIIndex]);

            const labelList& eJCells(mesh_.edgeCells()[edgeJIndex]);

            forAll(eICells, eIC)
            {
                const label cellI(eICells[eIC]);

                forAll(eJCells, eJC)
                {
                    const label cellJ(eJCells[eJC]);

                    if (cellJ > cellI)
                    {
                        if
                        (
                            findIndex(directInteractionList[cellI], cellJ)
                         == -1
                        )
                        {
                            directInteractionList[cellI].append(cellJ);
                        }
                    }

                    if (cellI > cellJ)
                    {
                        if
                        (
                            findIndex(directInteractionList[cellJ], cellI)
                         == -1
                        )
                        {
                            directInteractionList[cellJ].append(cellI);
                        }
                    }
                }
            }
        }
    }
}


// label pointJIndex;
//
// forAll (mesh_.points(), pointIIndex)
// {
//     const point& ptI
//     (
//         mesh_.points()[pointIIndex]
//     );
//
//     for
//     (
//         pointJIndex = pointIIndex;
//         pointJIndex != mesh_.points().size();
//         ++pointJIndex
//     )
//     {
//         const point& ptJ
//         (
//             mesh_.points()[pointJIndex]
//         );
//
//         if (magSqr(ptI - ptJ) <= rCutMaxSqr)
//         {
//             const labelList& ptICells
//             (
//                 mesh_.pointCells()[pointIIndex]
//             );
//
//             const labelList& ptJCells
//             (
//                 mesh_.pointCells()[pointJIndex]
//             );
//
//             forAll(ptICells, pIC)
//             {
//                 const label cellI(ptICells[pIC]);
//
//                 forAll(ptJCells, pJC)
//                 {
//                     const label cellJ(ptJCells[pJC]);
//
//                     if (cellJ > cellI)
//                     {
//                         if(findIndex(directInteractionList[cellI], cellJ) == -1)
//                         {
//                             directInteractionList[cellI].append(cellJ);
//                         }
//                     }
//
//                     if (cellI > cellJ)
//                     {
//                         if(findIndex(directInteractionList[cellJ], cellI) == -1)
//                         {
//                             directInteractionList[cellJ].append(cellI);
//                         }
//                     }
//                 }
//             }
//         }
//     }
// }
