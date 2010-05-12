labelVector iN(0,0,0);

vector gap = (vector::one)*pow((numberDensity/2.0),-(1.0/3.0));

#include "origin.H"

// Info<< "gap = " << gap << endl;

// Special treatment is required for the first position, i.e. iteration zero.

if (n == 0)
{
    latticePosition.x() = (iN.x()*gap.x());

    latticePosition.y() = (iN.y()*gap.y());

    latticePosition.z() = (iN.z()*gap.z());

    // Placing 2 molecules in each unit cell, using the algorithm from
    // D. Rapaport, The Art of Molecular Dynamics Simulation, 2nd Ed, p68

    for (label iU = 0; iU < 2; iU++)
    {
        vector unitCellLatticePosition = latticePosition;

        if (iU == 1)
        {
            unitCellLatticePosition += 0.5 * gap;
        }

        if (originSpecifies == "corner")
        {
            unitCellLatticePosition -= 0.25*gap;
        }

//         Info << nl << n << ", " << unitCellLatticePosition;

        globalPosition =
            origin + transform(latticeToGlobal,unitCellLatticePosition);

        partOfLayerInBounds = mesh_.bounds().contains(globalPosition);

        if
        (
             findIndex(mesh_.cellZones()[cZ], mesh_.findCell(globalPosition))
           != -1
        )
        {
            molsPlacedThisIteration++;

            initialPositions.append(globalPosition);

            initialCelli.append(mesh_.findCell(globalPosition));
        }
    }
}
else
{
    // Place top and bottom caps.

    for (iN.z() = -n; iN.z() <= n; iN.z() += 2*n)
    {
        for (iN.y() = -n; iN.y() <= n; iN.y()++)
        {
            for (iN.x() = -n; iN.x() <= n; iN.x()++)
            {
                latticePosition.x() = (iN.x() * gap.x());

                latticePosition.y() = (iN.y() * gap.y());

                latticePosition.z() = (iN.z() * gap.z());

                for (label iU = 0; iU < 2; iU++)
                {
                    vector unitCellLatticePosition = latticePosition;

                    if (iU == 1)
                    {
                        unitCellLatticePosition += 0.5*gap;
                    }

                    if(originSpecifies == "corner")
                    {
                        unitCellLatticePosition -= 0.25*gap;
                    }

//                     Info << nl << iN << ", " << unitCellLatticePosition;

                    globalPosition =
                        origin
                      + transform(latticeToGlobal,unitCellLatticePosition);

                    partOfLayerInBounds =
                        mesh_.bounds().contains(globalPosition);

                    if
                    (
                        findIndex
                        (
                             mesh_.cellZones()[cZ],
                             mesh_.findCell(globalPosition)
                        )
                     != -1)
                    {
                        molsPlacedThisIteration++;

                        initialPositions.append(globalPosition);

                        initialCelli.append(mesh_.findCell(globalPosition));
                    }
                }
            }
        }
    }

//     Placing sides

    for (iN.z() = -(n-1); iN.z() <= (n-1); iN.z()++)
    {
        for (label iR = 0; iR <= 2*n -1; iR++)
        {
            latticePosition.x() = (n*gap.x());

            latticePosition.y() = ((-n + (iR + 1))*gap.y());

            latticePosition.z() = (iN.z() * gap.z());

            for (label iK = 0; iK < 4; iK++)
            {
                for (label iU = 0; iU < 2; iU++)
                {
                    vector unitCellLatticePosition = latticePosition;

                    if (iU == 1)
                    {
                        unitCellLatticePosition += 0.5 * gap;
                    }

                    if (originSpecifies == "corner")
                    {
                        unitCellLatticePosition -= 0.25*gap;
                    }

                    globalPosition =
                        origin
                      + transform(latticeToGlobal,unitCellLatticePosition);

                    partOfLayerInBounds =
                        mesh_.bounds().contains(globalPosition);

                    if
                    (
                        findIndex
                         (
                             mesh_.cellZones()[cZ],
                             mesh_.findCell(globalPosition)
                         )
                      != -1
                    )
                    {
                        molsPlacedThisIteration++;

                        initialPositions.append(globalPosition);

                        initialCelli.append(mesh_.findCell(globalPosition));
                    }
                }

                latticePosition =
                    vector
                    (
                      - latticePosition.y(),
                        latticePosition.x(),
                        latticePosition.z()
                    );
            }
        }
    }
}
