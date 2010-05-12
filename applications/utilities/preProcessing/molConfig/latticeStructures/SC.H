labelVector iN(0,0,0);

vector gap = (vector::one)*pow(numberDensity, -(1.0/3.0));

#include "origin.H"

// Info<< "gap = " << gap << endl;

// Special treatment is required for the first position, i.e. iteration zero.

if (n == 0)
{
    latticePosition = vector::zero;

    if (originSpecifies == "corner")
    {
        latticePosition += 0.5*gap;
    }

    globalPosition = origin + transform(latticeToGlobal,latticePosition);

    partOfLayerInBounds = mesh_.bounds().contains(globalPosition);

    if (findIndex(mesh_.cellZones()[cZ], mesh_.findCell(globalPosition)) != -1)
    {
        molsPlacedThisIteration++;

        initialPositions.append(globalPosition);

        initialCelli.append(mesh_.findCell(globalPosition));
    }
}
else
{
    for (iN.z() = -n; iN.z() <= n; iN.z() += 2*n)
    {
        for (iN.y() = -n; iN.y() <= n; iN.y()++)
        {
            for (iN.x() = -n; iN.x() <= n; iN.x()++)
            {
                latticePosition.x() = (iN.x() * gap.x());

                latticePosition.y() = (iN.y() * gap.y());

                latticePosition.z() = (iN.z() * gap.z());

                if (originSpecifies == "corner")
                {
                    latticePosition += 0.5*gap;
                }

                globalPosition =
                    origin + transform(latticeToGlobal,latticePosition);

                partOfLayerInBounds = mesh_.bounds().contains(globalPosition);

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
        }
    }

    tensor quarterRotate(EulerCoordinateRotation(-90, 0, 0, true).R());

    iN.x() = n;
    for (iN.z() = -(n-1); iN.z() <= (n-1); iN.z()++)
    {
        for (iN.y() = -(n-1); iN.y() <= n; iN.y()++)
        {
            latticePosition.x() = (iN.x()*gap.x());

            latticePosition.y() = (iN.y()*gap.y());

            latticePosition.z() = (iN.z()*gap.z());

            for (label iR = 0; iR < 4; iR++)
            {
                vector offsetCorrectedLatticePosition = latticePosition;

                if (originSpecifies == "corner")
                {
                    offsetCorrectedLatticePosition += 0.5*gap;
                }

                globalPosition =
                    origin
                  + transform(latticeToGlobal,offsetCorrectedLatticePosition);

                partOfLayerInBounds = mesh_.bounds().contains(globalPosition);

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

                latticePosition = transform(quarterRotate,latticePosition);
            }
        }
    }
}

