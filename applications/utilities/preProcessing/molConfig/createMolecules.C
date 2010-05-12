/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "molConfig.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::molConfig::createMolecules()
{
    Info<< nl << "Creating molecules from zone specifications\n" << endl;

    DynamicList<vector> initialPositions(0);

    DynamicList<label> initialIds(0);

    DynamicList<scalar> initialMasses(0);

    DynamicList<label> initialCelli(0);

    DynamicList<vector> initialVelocities(0);

    DynamicList<vector> initialAccelerations(0);

    DynamicList<label> initialTethered(0);

    DynamicList<vector> initialTetherPositions(0);

    label totalMols = 0;

    label idAssign;

    Random rand(clock::getTime());

// * * * * * * * * Building the IdList * * * * * * * * * //

//Notes: - each processor will have an identical idList_.
//       - The order of id's inside the idList_ depends on the order
//         of subDicts inside the molConigDict.

    Info<< "Building the idList: " ;

    forAll(molConfigDescription_.toc(), cZs)
    {
        word subDictName (molConfigDescription_.toc()[cZs]);

        word iD (molConfigDescription_.subDict(subDictName).lookup("id"));

        if (findIndex(idList_,iD) == -1)
        {
            idList_.append(iD);
        }
    }

    forAll(idList_, i)
    {
        Info << " " << idList_[i];
    }

    Info << nl << endl;

// * * * * * * * * Filling the Mesh * * * * * * * * * //

    const cellZoneMesh& cellZoneI = mesh_.cellZones();

    if (cellZoneI.size())
    {
        Info<< "Filling the zones with molecules." << nl << endl;
    }
    else
    {
        FatalErrorIn("void createMolecules()\n")
            << "No cellZones found in mesh description."
            << abort(FatalError);
    }

    forAll (cellZoneI, cZ)
    {
        if (cellZoneI[cZ].size())
        {
            if (!molConfigDescription_.found(cellZoneI[cZ].name()))
            {
                Info << "Zone specification subDictionary: "
                    << cellZoneI[cZ].name() << " not found." << endl;
            }
            else
            {
                label n = 0;

                label totalZoneMols = 0;

                // Garbage: uninitialised data.  HJ, 27/Nov/2009
                label molsPlacedThisIteration;

#               include "readZoneSubDict.H"

                idAssign = findIndex(idList_,id);

#               include "startingPoint.H"

                // Continue trying to place molecules as long as at
                // least one molecule is placed in each iteration.
                // The "|| totalZoneMols == 0" condition means that the
                // algorithm will continue if the origin is outside the
                // zone - it will cause an infinite loop if no molecules
                // are ever placed by the algorithm.

                if (latticeStructure != "empty")
                {
                    while
                    (
                        molsPlacedThisIteration != 0
                     || totalZoneMols == 0
                    )
                    {
                        molsPlacedThisIteration = 0;

                        bool partOfLayerInBounds = false;

#                       include "createPositions.H"

                        if
                        (
                            totalZoneMols == 0
                         && !partOfLayerInBounds
                        )
                        {
                            WarningIn("molConfig::createMolecules()")
                                << "A whole layer of unit cells was placed "
                                << "outside the bounds of the mesh, but no "
                                << "molecules have been placed in zone '"
                                << cellZoneI[cZ].name()
                                << "'.  This is likely to be because the zone "
                                << "has few cells ("
                                << cellZoneI[cZ].size()
                                << " in this case) and no lattice position "
                                << "fell inside them.  "
                                << "Aborting filling this zone."
                                << endl;

                            break;
                        }

                        totalZoneMols += molsPlacedThisIteration;

                        n++;
                    }

                    label molN;

                    for
                    (
                        molN = totalMols;
                        molN < totalMols + totalZoneMols;
                        molN++
                    )
                    {
                        initialIds.append(idAssign);

                        initialMasses.append(mass);

                        initialAccelerations.append(vector::zero);

                        if (tethered)
                        {
                            initialTethered.append(1);

                            initialTetherPositions.append
                            (
                                initialPositions[molN]
                            );
                        }

                        else
                        {
                            initialTethered.append(0);

                            initialTetherPositions.append(vector::zero);
                        }
                    }

#                   include "createVelocities.H"

#                   include "correctVelocities.H"

                }

                totalMols += totalZoneMols;
            }
        }
    }

    idList_.shrink();

    positions_ = initialPositions;

    positions_.setSize(initialPositions.size());

    id_ = initialIds;

    id_.setSize(initialIds.size());

    mass_ = initialMasses;

    mass_.setSize(initialMasses.size());

    cells_ = initialCelli;

    cells_.setSize(initialCelli.size());

    U_ = initialVelocities;

    U_.setSize(initialVelocities.size());

    A_ = initialAccelerations;

    A_.setSize(initialAccelerations.size());

    tethered_ = initialTethered;

    tethered_.setSize(initialTethered.size());

    tetherPositions_ = initialTetherPositions;

    tetherPositions_.setSize(initialTetherPositions.size());

    nMol_ = totalMols;
}


// ************************************************************************* //
