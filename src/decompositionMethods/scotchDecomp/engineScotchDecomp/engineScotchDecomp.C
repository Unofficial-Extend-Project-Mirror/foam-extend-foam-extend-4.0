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

#include "engineScotchDecomp.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(engineScotchDecomp, 0);

    addToRunTimeSelectionTable
    (
        scotchDecomp,
        engineScotchDecomp,
        dictionaryMesh
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::engineScotchDecomp::engineScotchDecomp
(
    const dictionary& decompositionDict,
    const polyMesh& mesh
)
:
    scotchDecomp(decompositionDict, mesh),
    dict_(decompositionDict.subDict(typeName + "Coeffs")),
    slidingPatchPairs_(dict_.lookup("slidingPatchPairs")),
    expandSliding_(dict_.lookup("expandSliding"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::engineScotchDecomp::decompose
(
    const pointField& points,
    const scalarField& pointWeights
)
{
    if (points.size() != mesh().nCells())
    {
        FatalErrorIn
        (
            "engineScotchDecomp::decompose\n"
            "(\n"
            "    const pointField&,\n"
            "    const scalarField&\n"
            ")"
        )   << "Can use this decomposition method only for the whole mesh"
            << endl
            << "and supply one coordinate (cellCentre) for every cell." << endl
            << "The number of coordinates " << points.size() << endl
            << "The number of cells in the mesh " << mesh().nCells()
            << exit(FatalError);
    }

    // Create clustering to coarse level

    labelList fineToCoarse(mesh().nCells(), -1);

    // Mask all cells in the cylinder layering region with true
    // Used in two-pass decomposition later
    boolList cylinderMask(mesh().nCells(), false);
    label nClusters = 0;

    // Locate piston patch and cluster up colums, using opposite faces
    label pistonPatchID = mesh().boundaryMesh().findPatchID("piston");

    // Go through the sliding pairs and mark clustering
    forAll (slidingPatchPairs_, pairI)
    {
        // Locate patch and cluster up colums, using opposite faces
        label firstPatchID = mesh().boundaryMesh().findPatchID
        (
            slidingPatchPairs_[pairI].first()
        );

        label secondPatchID = mesh().boundaryMesh().findPatchID
        (
            slidingPatchPairs_[pairI].second()
        );

        if (firstPatchID == -1 || secondPatchID == -1)
        {
            FatalErrorIn
            (
                "labelList engineScotchDecomp::decompose\n"
                "(\n"
                "    const pointField& points,\n"
                "    const scalarField& pointWeights\n"
                ")"
            )   << "Cannot find sliding patch pair "
                << slidingPatchPairs_[pairI]
                << abort(FatalError);
        }

        // Put all cells next to the patch into a cluster
        if (expandSliding_)
        {
            // Use point-cell addressing from patch points
            const labelListList& pc = mesh().pointCells();

            // First side
            const labelList& mp1 =
                mesh().boundaryMesh()[firstPatchID].meshPoints();

            forAll (mp1, pointI)
            {
                const labelList& curCells = pc[mp1[pointI]];

                forAll (curCells, cellI)
                {
                    fineToCoarse[curCells[cellI]] = nClusters;
                    cylinderMask[curCells[cellI]] = true;
                }
            }

            // Second side
            {
                const labelList& mp2 =
                    mesh().boundaryMesh()[secondPatchID].meshPoints();

                forAll (mp2, pointI)
                {
                    const labelList& curCells = pc[mp2[pointI]];

                    forAll (curCells, cellI)
                    {
                        fineToCoarse[curCells[cellI]] = nClusters;
                        cylinderMask[curCells[cellI]] = true;
                    }
                }
            }
        }
        else
        {
            // First side
            {
                const labelList& fc1 =
                    mesh().boundaryMesh()[firstPatchID].faceCells();

                forAll (fc1, fcI)
                {
                    fineToCoarse[fc1[fcI]] = nClusters;
                    cylinderMask[fc1[fcI]] = true;
                }
            }

            // Second side
            {
                const labelList& fc2 =
                    mesh().boundaryMesh()[secondPatchID].faceCells();

                forAll (fc2, fcI)
                {
                    fineToCoarse[fc2[fcI]] = nClusters;
                    cylinderMask[fc2[fcI]] = true;
                }
            }
        }

        nClusters++;
    }

    if (pistonPatchID > -1)
    {
        // Found piston patch
        const faceList& f = mesh().allFaces();
        const cellList& c = mesh().cells();
        const labelList& owner = mesh().faceOwner();
        const labelList& neighbour = mesh().faceNeighbour();

        const labelList& faceCells =
            mesh().boundaryMesh()[pistonPatchID].faceCells();

        forAll (faceCells, faceI)
        {
            // Get face index
            label curFaceNo = faceI
                + mesh().boundaryMesh()[pistonPatchID].start();

            // Get cell index
            label curCellNo = faceCells[faceI];

            // Mark cell to cluster
            if (fineToCoarse[curCellNo] < 0)
            {
                // New cluster
                fineToCoarse[curCellNo] = nClusters;
                cylinderMask[curCellNo] = true;

                for(;;)
                {
                    // Attempt to find the next face and cell
                    curFaceNo = c[curCellNo].opposingFaceLabel(curFaceNo, f);

                    if (curFaceNo > -1)
                    {
                        // Face found, try for a cell
                        if (curFaceNo < mesh().nInternalFaces())
                        {
                            if (owner[curFaceNo] == curCellNo)
                            {
                                curCellNo = neighbour[curFaceNo];
                            }
                            else if (neighbour[curFaceNo] == curCellNo)
                            {
                                curCellNo = owner[curFaceNo];
                            }
                            else
                            {
                                // Error in layering.  Should never happen
                                break;
                            }

                            // Mark cell to cluster
                            fineToCoarse[curCellNo] = nClusters;
                            cylinderMask[curCellNo] = true;
                        }
                        else
                        {
                            // Hit boundary face
                            break;
                        }
                    }
                    else
                    {
                        // Cannot find opposing face: out of prismatic region
                        break;
                    }
                }

                // Next cluster
                nClusters++;
            }
        }
    }

    // Count cylinder cells from mask
    label nCylinderCells = 0;

    forAll (cylinderMask, cellI)
    {
        if (cylinderMask[cellI]) nCylinderCells++;
    }

    label nStaticCells = mesh().nCells() - nCylinderCells;
    label nCylClusters = nClusters;

    // Visit all unmarked cells and put them into single clusters
    forAll (fineToCoarse, cellI)
    {
        if (fineToCoarse[cellI] == -1)
        {
            fineToCoarse[cellI] = nClusters;
            nClusters++;
        }
    }

    label nStaticClusters = nClusters - nCylClusters;

    Info<< "Number of cells: " << mesh().nCells()
        << ", in cylinder + sliding: " << nCylinderCells
        << ", in static part: " << nStaticCells << nl
        << "Number of cylinder clusters " << nCylClusters
        << ", static clusters " << nStaticClusters
        << ", total clusters " << nClusters << endl;

    // Mark-up complete.  Create point centres and weights for all clusters
    vectorField clusterCentres(nClusters, vector::zero);

    // Stabilise cluster volumes in case a cluster ends up empty
    // It is possible to have empty clusters without connectivity
    scalarField clusterVols(nClusters, SMALL);
    scalarField clusterWeights(nClusters, 0);

    const vectorField& centres = mesh().cellCentres();
    const scalarField& vols = mesh().cellVolumes();

    forAll (fineToCoarse, cellI)
    {
        const label& curCoarse = fineToCoarse[cellI];

        clusterCentres[curCoarse] += centres[cellI]*vols[cellI];
        clusterVols[curCoarse] += vols[cellI];
        clusterWeights[curCoarse] += 1;
    }

    clusterCentres /= clusterVols;

    // Execute decomposition, first on cylinder layering zone, then on the rest

    // Collect cell-cells in cylinder layering zone and the rest

    // Create and cellCells hash lookup on two pieces

    // Note: cell clusters come first and they will be done
    // on a shortened list.  Static clusters need to be renumbered by
    // throwing away the first part of the list

    List<labelHashSet> cylCellCellsHash(nCylClusters);
    List<labelHashSet> staticCellCellsHash(nStaticClusters);

    const labelListList cellCells = mesh().cellCells();

    forAll (cellCells, cellI)
    {
        const labelList& curCC = cellCells[cellI];

        label curCluster = fineToCoarse[cellI];

        if (cylinderMask[cellI])
        {
            labelHashSet& curCylAddr = cylCellCellsHash[curCluster];

            // Collect neighbour cluster addressing
            forAll (curCC, neiI)
            {
                // Add neighbour if marked
                if (cylinderMask[curCC[neiI]])
                {
                    label nbrCluster = fineToCoarse[curCC[neiI]];

                    if (nbrCluster != curCluster)
                    {
                        if (!curCylAddr.found(nbrCluster))
                        {
                            curCylAddr.insert(nbrCluster);
                        }
                    }
                }
            }
        }
        else
        {
            // Offset index
            curCluster -= nCylClusters;

            labelHashSet& curStaticAddr = staticCellCellsHash[curCluster];

            forAll (curCC, neiI)
            {
                // Add neighbour if marked
                if (!cylinderMask[curCC[neiI]])
                {
                    label nbrCluster = fineToCoarse[curCC[neiI]]
                        - nCylClusters;

                    if (nbrCluster != curCluster)
                    {
                        if (!curStaticAddr.found(nbrCluster))
                        {
                            curStaticAddr.insert(nbrCluster);
                        }
                    }
                }
            }
        }
    }

    // Pack cellCells on the cylinder
    labelListList cylCellCells(nCylClusters);

    forAll (cylCellCellsHash, clusterI)
    {
        cylCellCells[clusterI] = cylCellCellsHash[clusterI].toc();
    }

    // Decompose cylinder: size of list equals the number of clusters
    // in the cylinder region
    vectorField clusterCentresCyl = clusterCentres;
    scalarField clusterWeightsCyl = clusterWeights;

    clusterCentresCyl.setSize(nCylClusters);
    clusterWeightsCyl.setSize(nCylClusters);

    labelList cylDecomp = scotchDecomp::decompose
    (
        cylCellCells,
        clusterCentresCyl,
        clusterWeightsCyl
    );

    // Decompose static: size of list equals the number of clusters

    labelListList staticCellCells(nStaticClusters);

    forAll (staticCellCellsHash, clusterI)
    {
        staticCellCells[clusterI] = staticCellCellsHash[clusterI].toc();
    }

    vectorField clusterCentresStatic(nStaticClusters);
    scalarField clusterWeightsStatic(nStaticClusters);

    forAll (clusterCentresStatic, i)
    {
        clusterCentresStatic[i] = clusterCentres[nCylClusters + i];
        clusterWeightsStatic[i] = clusterWeights[nCylClusters + i];
    }

    labelList staticDecomp = scotchDecomp::decompose
    (
        staticCellCells,
        clusterCentresStatic,
        clusterWeightsStatic
    );

    // Reconstruct final decomposition

    labelList finalDecomp(mesh().nCells(), -1);

    forAll (cylinderMask, cellI)
    {
        if (cylinderMask[cellI])
        {
            // Cylinder cell
            finalDecomp[cellI] = cylDecomp[fineToCoarse[cellI]];
        }
        else
        {
            // Static cell
            finalDecomp[cellI] =
                staticDecomp[fineToCoarse[cellI] - nCylClusters];
        }
    }

    return finalDecomp;
}


// ************************************************************************* //
