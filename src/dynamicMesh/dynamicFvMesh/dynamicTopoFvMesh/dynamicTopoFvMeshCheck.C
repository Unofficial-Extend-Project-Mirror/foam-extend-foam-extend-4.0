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

Class
    dynamicTopoFvMesh

Description
    Functions specific to connectivity checking and debugging

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "dynamicTopoFvMesh.H"

#include "IOmanip.H"
#include "volFields.H"
#include "triPointRef.H"
#include "tetPointRef.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Compute mesh-quality, and return true if no slivers are present
bool dynamicTopoFvMesh::meshQuality
(
    bool outputOption
)
{
    // Compute statistics on the fly
    label nCells = 0, minCell = -1;
    scalar maxQuality = -GREAT;
    scalar minQuality =  GREAT;
    scalar cQuality, meanQuality = 0.0;

    // Track slivers
    bool sliversAbsent = true;
    thresholdSlivers_.clear();

    // Loop through all cells in the mesh and compute cell quality
    forAll(cells_, cellI)
    {
        const cell& cellToCheck = cells_[cellI];

        if (cellToCheck.empty())
        {
            continue;
        }

        if (twoDMesh_)
        {
            // Assume XY plane here
            vector n = vector(0,0,1);

            // Get a triangular boundary face
            forAll(cellToCheck, faceI)
            {
                const face& faceToCheck = faces_[cellToCheck[faceI]];

                if (faceToCheck.size() == 3)
                {
                    triPointRef tpr
                    (
                        points_[faceToCheck[0]],
                        points_[faceToCheck[1]],
                        points_[faceToCheck[2]]
                    );

                    // Assume centre-plane passes through origin
                    cQuality =
                    (
                        tpr.quality() *
                        (
                            Foam::sign
                            (
                                tpr.normal() &
                                ((tpr.centre() & n) * n)
                            )
                        )
                    );

                    break;
                }
            }
        }
        else
        {
            const label bfIndex = cellToCheck[0];
            const label cfIndex = cellToCheck[1];

            const face& baseFace = faces_[bfIndex];
            const face& checkFace = faces_[cfIndex];

            // Get the fourth point
            label apexPoint =
            (
                meshOps::findIsolatedPoint(baseFace, checkFace)
            );

            // Compute cell quality
            if (owner_[bfIndex] == cellI)
            {
                cQuality =
                (
                    tetMetric_
                    (
                        points_[baseFace[2]],
                        points_[baseFace[1]],
                        points_[baseFace[0]],
                        points_[apexPoint]
                    )
                );
            }
            else
            {
                cQuality =
                (
                    tetMetric_
                    (
                        points_[baseFace[0]],
                        points_[baseFace[1]],
                        points_[baseFace[2]],
                        points_[apexPoint]
                    )
                );
            }
        }

        // Update statistics
        maxQuality = Foam::max(cQuality, maxQuality);

        if (cQuality < minQuality)
        {
            minQuality = cQuality;
            minCell = cellI;
        }

        meanQuality += cQuality;
        nCells++;

        // Add to the list of slivers
        if ((cQuality < sliverThreshold_) && (cQuality > 0.0))
        {
            thresholdSlivers_.insert(cellI, cQuality);
        }
    }

    if (thresholdSlivers_.size())
    {
        sliversAbsent = false;
    }

    // Reduce across processors.
    reduce(sliversAbsent, andOp<bool>());

    // Output statistics:
    if (outputOption || (debug > 0))
    {
        // Reduce statistics across processors.
        reduce(minQuality, minOp<scalar>());
        reduce(maxQuality, maxOp<scalar>());
        reduce(meanQuality, sumOp<scalar>());
        reduce(nCells, sumOp<label>());

        if (minQuality < 0.0)
        {
            WarningIn
            (
                "bool dynamicTopoFvMesh::meshQuality"
                "(bool outputOption)"
            )
                << nl
                << "Minimum cell quality is: " << minQuality
                << " at cell: " << minCell
                << endl;
        }

        Info << endl;
        Info << " ~~~ Mesh Quality Statistics ~~~ " << endl;
        Info << " Min: " << minQuality << endl;
        Info << " Max: " << maxQuality << endl;
        Info << " Mean: " << meanQuality/nCells << endl;
        Info << " Cells: " << nCells << endl;
        Info << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << nl << endl;
    }

    return sliversAbsent;
}


// Utility to check whether points of an edge lie on a boundary.
const FixedList<bool,2>
dynamicTopoFvMesh::checkEdgeBoundary
(
    const label eIndex
) const
{
    FixedList<bool,2> edgeBoundary(false);

    const edge& edgeToCheck = edges_[eIndex];

    // Loop through edges connected to both points,
    // and check if any of them lie on boundaries.
    // Used to ensure that collapses happen towards boundaries.
    forAll(edgeToCheck, pointI)
    {
        const labelList& pEdges = pointEdges_[edgeToCheck[pointI]];

        forAll(pEdges, edgeI)
        {
            // Determine the patch this edge belongs to
            if (whichEdgePatch(pEdges[edgeI]) > -1)
            {
                edgeBoundary[pointI] = true;
                break;
            }
        }
    }

    return edgeBoundary;
}


// Check whether the given edge is on a bounding curve
bool dynamicTopoFvMesh::checkBoundingCurve(const label eIndex) const
{
    // Internal edges don't count
    label edgePatch = -1;

    if ((edgePatch = whichEdgePatch(eIndex)) < 0)
    {
        return false;
    }
    else
    {
        // Check whether this edge shouldn't be swapped
        if (findIndex(noSwapPatchIDs_, edgePatch) > -1)
        {
            return true;
        }
    }

    // Check if two boundary faces lie on different face-patches
    FixedList<vector, 2> fNorm(vector::zero);
    label fPatch, firstPatch = -1, secondPatch = -1, count = 0;
    const labelList& edgeFaces = edgeFaces_[eIndex];

    forAll(edgeFaces, faceI)
    {
        if ((fPatch = whichPatch(edgeFaces[faceI])) > -1)
        {
            // Obtain the normal.
            fNorm[count] = faces_[edgeFaces[faceI]].normal(points_);

            // Normalize it.
            fNorm[count] /= mag(fNorm[count]) + VSMALL;

            count++;

            if (firstPatch == -1)
            {
                firstPatch = fPatch;
            }
            else
            {
                secondPatch = fPatch;
                break;
            }
        }
    }

    scalar deviation = (fNorm[0] & fNorm[1]);

    // Check if the swap-curvature is too high
    if (mag(deviation) < swapDeviation_)
    {
        return true;
    }

    // Check if the edge borders two different patches
    if (firstPatch != secondPatch)
    {
        return true;
    }

    // Not on a bounding curve
    return false;
}


// Check triangulation quality for an edge index
bool dynamicTopoFvMesh::checkQuality
(
    const label eIndex,
    const labelList& m,
    const PtrList<scalarListList>& Q,
    const scalar minQuality,
    const label checkIndex
) const
{
    bool myResult = false;

    // Non-coupled check
    if (Q[checkIndex][0][m[checkIndex]-1] > minQuality)
    {
        myResult = true;

        if (debug > 2)
        {
            Info << nl << nl
                 << " eIndex: " << eIndex
                 << " minQuality: " << minQuality
                 << " newQuality: " << Q[checkIndex][0][m[checkIndex]-1]
                 << endl;
        }
    }

    return myResult;
}


// Print out tables for debugging
void dynamicTopoFvMesh::printTables
(
    const labelList& m,
    const PtrList<scalarListList>& Q,
    const PtrList<labelListList>& K,
    const label checkIndex
) const
{
    Info << "m: " << m[checkIndex] << endl;

    // Print out Q
    Info << "===" << endl;
    Info << " Q " << endl;
    Info << "===" << endl;

    Info << "   ";

    for (label j = 0; j < m[checkIndex]; j++)
    {
        Info << setw(12) << j;
    }

    Info << nl;

    for (label j = 0; j < m[checkIndex]; j++)
    {
        Info << "-------------";
    }

    Info << nl;

    for (label i = 0; i < (m[checkIndex]-2); i++)
    {
        Info << i << ": ";

        for (label j = 0; j < m[checkIndex]; j++)
        {
            Info << setw(12) << Q[checkIndex][i][j];
        }

        Info << nl;
    }

    // Print out K
    Info << "===" << endl;
    Info << " K " << endl;
    Info << "===" << endl;

    Info << "   ";

    for (label j = 0; j < m[checkIndex]; j++)
    {
        Info << setw(12) << j;
    }

    Info << nl;

    for (label i = 0; i < (m[checkIndex]-2); i++)
    {
        Info << i << ": ";

        for (label j = 0; j < m[checkIndex]; j++)
        {
            Info << setw(12) << K[checkIndex][i][j];
        }

        Info << nl;
    }

    Info << endl;
}


// Check old-volumes for an input triangulation
bool dynamicTopoFvMesh::checkTriangulationVolumes
(
    const label eIndex,
    const labelList& hullVertices,
    const labelListList& triangulations
) const
{
    label m = hullVertices.size();
    scalar tetVol = 0.0;

    const edge& edgeToCheck = edges_[eIndex];

    for (label i = 0; i < (m-2); i++)
    {
        // Compute volume for the upper-half
        tetVol =
        (
            tetPointRef
            (
                oldPoints_[hullVertices[triangulations[0][i]]],
                oldPoints_[hullVertices[triangulations[1][i]]],
                oldPoints_[hullVertices[triangulations[2][i]]],
                oldPoints_[edgeToCheck[0]]
            ).mag()
        );

        if (tetVol < 0.0)
        {
            if (debug > 2)
            {
                InfoIn
                (
                    "bool dynamicTopoFvMesh::checkTriangulationVolumes\n"
                    "(\n"
                    "    const label eIndex,\n"
                    "    const labelList& hullVertices,\n"
                    "    const labelListList& triangulations\n"
                    ") const\n"
                )
                    << "Swap sequence leads to negative old-volumes." << nl
                    << "Edge: " << edgeToCheck << nl
                    << "using Points: " << nl
                    << oldPoints_[hullVertices[triangulations[0][i]]] << nl
                    << oldPoints_[hullVertices[triangulations[1][i]]] << nl
                    << oldPoints_[hullVertices[triangulations[2][i]]] << nl
                    << oldPoints_[edgeToCheck[0]] << endl;
            }

            return true;
        }

        tetVol =
        (
            tetPointRef
            (
                oldPoints_[hullVertices[triangulations[2][i]]],
                oldPoints_[hullVertices[triangulations[1][i]]],
                oldPoints_[hullVertices[triangulations[0][i]]],
                oldPoints_[edgeToCheck[1]]
            ).mag()
        );

        if (tetVol < 0.0)
        {
            if (debug > 2)
            {
                InfoIn
                (
                    "bool dynamicTopoFvMesh::checkTriangulationVolumes\n"
                    "(\n"
                    "    const label eIndex,\n"
                    "    const labelList& hullVertices,\n"
                    "    const labelListList& triangulations\n"
                    ") const\n"
                )
                    << "Swap sequence leads to negative old-volumes." << nl
                    << "Edge: " << edgeToCheck << nl
                    << "using Points: " << nl
                    << oldPoints_[hullVertices[triangulations[2][i]]] << nl
                    << oldPoints_[hullVertices[triangulations[1][i]]] << nl
                    << oldPoints_[hullVertices[triangulations[0][i]]] << nl
                    << oldPoints_[edgeToCheck[1]] << endl;
            }

            return true;
        }
    }

    return false;
}


// Output an entity as a VTK file
void dynamicTopoFvMesh::writeVTK
(
    const word& name,
    const label entity,
    const label primitiveType,
    const bool useOldConnectivity,
    const bool useOldPoints
) const
{
    writeVTK
    (
        name,
        labelList(1, entity),
        primitiveType,
        useOldConnectivity,
        useOldPoints
    );
}


// Output a list of primitives as a VTK file.
//  - primitiveType is:
//      0: List of points
//      1: List of edges
//      2: List of faces
//      3: List of cells
void dynamicTopoFvMesh::writeVTK
(
    const word& name,
    const labelList& cList,
    const label primitiveType,
    const bool useOldConnectivity,
    const bool useOldPoints
) const
{
    if (useOldPoints)
    {
        if (useOldConnectivity)
        {
            // Use points from polyMesh
            meshOps::writeVTK
            (
                (*this),
                name,
                cList,
                primitiveType,
                polyMesh::points(),
                polyMesh::edges(),
                polyMesh::faces(),
                polyMesh::cells(),
                polyMesh::faceOwner()
            );
        }
        else
        {
            meshOps::writeVTK
            (
                (*this),
                name,
                cList,
                primitiveType,
                oldPoints_,
                edges_,
                faces_,
                cells_,
                owner_
            );
        }
    }
    else
    {
        meshOps::writeVTK
        (
            (*this),
            name,
            cList,
            primitiveType,
            points_,
            edges_,
            faces_,
            cells_,
            owner_
        );
    }
}


// Return the status report interval
scalar dynamicTopoFvMesh::reportInterval() const
{
    // Default to 3 seconds
    scalar interval = 3.0;

    const dictionary& meshSubDict = dict_.subDict("dynamicTopoFvMesh");

    if (meshSubDict.found("reportInterval") || mandatory_)
    {
        interval = readScalar(meshSubDict.lookup("reportInterval"));

        // Prevent reports if necessary
        if (interval < VSMALL)
        {
            interval = GREAT;
        }
    }

    return interval;
}


// Check the state of connectivity lists
void dynamicTopoFvMesh::checkConnectivity(const label maxErrors) const
{
    label nFailedChecks = 0;

    messageStream ConnectivityWarning
    (
        "dynamicTopoFvMesh Connectivity Warning",
        messageStream::SERIOUS,
        maxErrors
    );

    // Check face-label ranges
    Info << "Checking index ranges...";

    forAll(edges_, edgeI)
    {
        const edge& curEdge = edges_[edgeI];

        if (curEdge == edge(-1, -1))
        {
            continue;
        }

        if
        (
            curEdge[0] < 0 || curEdge[0] > (points_.size()-1) ||
            curEdge[1] < 0 || curEdge[1] > (points_.size()-1)
        )
        {
            Pout << "Edge " << edgeI
                 << " contains vertex labels out of range: "
                 << curEdge
                 << " Max point index = " << (points_.size()-1) << endl;

            nFailedChecks++;

            ConnectivityWarning()
                 << nl << "Edge-point connectivity is inconsistent."
                 << endl;
        }

        // Check for unique point-labels
        if (curEdge[0] == curEdge[1])
        {
            Pout << "Edge " << edgeI
                 << " contains identical vertex labels: "
                 << curEdge
                 << endl;

            nFailedChecks++;

            ConnectivityWarning()
                 << nl << "Edge-point connectivity is inconsistent."
                 << endl;
        }
    }

    label allPoints = points_.size();
    labelList nPointFaces(allPoints, 0);

    forAll(faces_, faceI)
    {
        const face& curFace = faces_[faceI];

        if (curFace.empty())
        {
            continue;
        }

        if (min(curFace) < 0 || max(curFace) > (points_.size()-1))
        {
            Pout << "Face " << faceI
                 << " contains vertex labels out of range: "
                 << curFace
                 << " Max point index = " << (points_.size()-1) << endl;

            nFailedChecks++;

            ConnectivityWarning()
                 << nl << "Face-point connectivity is inconsistent."
                 << endl;
        }

        // Check for unique point-labels
        labelHashSet uniquePoints;

        forAll(curFace, pointI)
        {
            bool inserted = uniquePoints.insert(curFace[pointI]);

            if (!inserted)
            {
                Pout << "Face " << faceI
                     << " contains identical vertex labels: "
                     << curFace
                     << endl;

                nFailedChecks++;

                ConnectivityWarning()
                     << nl << "Face-point connectivity is inconsistent."
                     << endl;
            }
        }

        // Count faces per point
        forAll(curFace, pointI)
        {
            nPointFaces[curFace[pointI]]++;
        }

        // Ensure that cells on either side of this face
        // share just one face.
        if (neighbour_[faceI] > -1)
        {
            const cell& ownCell = cells_[owner_[faceI]];
            const cell& neiCell = cells_[neighbour_[faceI]];

            label nCommon = 0;

            forAll(ownCell, fi)
            {
                if (findIndex(neiCell, ownCell[fi]) > -1)
                {
                    nCommon++;
                }
            }

            if (nCommon != 1)
            {
                Pout << "Cells: " << nl
                     << '\t' << owner_[faceI] << ":: " << ownCell << nl
                     << '\t' << neighbour_[faceI] << " :: " << neiCell << nl
                     << " share multiple faces. "
                     << endl;

                nFailedChecks++;

                ConnectivityWarning()
                     << nl << "Cell-Face connectivity is inconsistent."
                     << endl;
            }
        }
    }

    label allFaces = faces_.size();
    labelList nCellsPerFace(allFaces, 0);

    forAll(cells_, cellI)
    {
        const cell& curCell = cells_[cellI];

        if (curCell.empty())
        {
            continue;
        }

        if (min(curCell) < 0 || max(curCell) > (faces_.size()-1))
        {
            Pout << "Cell " << cellI
                 << " contains vertex labels out of range: "
                 << curCell
                 << " Max point index = " << (faces_.size()-1) << endl;

            nFailedChecks++;

            ConnectivityWarning()
                 << nl << "Cell-Face connectivity is inconsistent."
                 << endl;
        }

        // Check for unique face-labels
        labelHashSet uniqueFaces;

        forAll(curCell, faceI)
        {
            bool inserted = uniqueFaces.insert(curCell[faceI]);

            if (!inserted)
            {
                Pout << "Cell " << cellI
                     << " contains identical face labels: "
                     << curCell
                     << endl;

                nFailedChecks++;

                ConnectivityWarning()
                     << nl << "Cell-Face connectivity is inconsistent."
                     << endl;
            }

            // Count cells per face
            nCellsPerFace[curCell[faceI]]++;
        }
    }

    Info << "Done." << endl;

    Info << "Checking face-cell connectivity...";

    forAll(nCellsPerFace, faceI)
    {
        if (nCellsPerFace[faceI] == 0)
        {
            // This might be a deleted face
            if (faceI < nOldFaces_)
            {
                if (reverseFaceMap_[faceI] == -1)
                {
                    continue;
                }
            }
            else
            {
                if (deletedFaces_.found(faceI))
                {
                    continue;
                }
            }

            // Looks like this is really an unused face.
            Pout << "Face " << faceI
                 << " :: " << faces_[faceI]
                 << " is unused. "
                 << endl;

            nFailedChecks++;

            ConnectivityWarning()
                 << nl << "Cell-Face connectivity is inconsistent."
                 << endl;
        }
        else
        if (nCellsPerFace[faceI] != 2 && whichPatch(faceI) == -1)
        {
            // Internal face is not shared by exactly two cells
            Pout << "Internal Face " << faceI
                 << " :: " << faces_[faceI]
                 << " is multiply connected." << nl
                 << " nCellsPerFace: " << nCellsPerFace[faceI]
                 << endl;

            nFailedChecks++;

            ConnectivityWarning()
                 << nl << "Cell-Face connectivity is inconsistent."
                 << endl;
        }
        else
        if (nCellsPerFace[faceI] != 1 && whichPatch(faceI) > -1)
        {
            // Boundary face is not shared by exactly one cell
            Pout << "Boundary Face " << faceI
                 << " :: " << faces_[faceI]
                 << " is multiply connected." << nl
                 << " nCellsPerFace: " << nCellsPerFace[faceI]
                 << endl;

            nFailedChecks++;

            ConnectivityWarning()
                 << nl << "Cell-Face connectivity is inconsistent."
                 << endl;
        }
    }

    Info << "Done." << endl;

    Info << "Checking for unused points...";

    forAll(nPointFaces, pointI)
    {
        if (nPointFaces[pointI] == 0)
        {
            // This might be a deleted point.
            if (pointI < nOldPoints_)
            {
                if (reversePointMap_[pointI] == -1)
                {
                    continue;
                }
            }
            else
            {
                if (deletedPoints_.found(pointI))
                {
                    continue;
                }
            }

            // Looks like this is really an unused point.
            Pout << "Point " << pointI << " is unused. " << endl;

            nFailedChecks++;

            ConnectivityWarning()
                 << nl << "Point-Face connectivity is inconsistent."
                 << endl;
        }
    }

    Info << "Done." << endl;

    Info << "Checking edge-face connectivity...";

    label allEdges = edges_.size();
    labelList nEdgeFaces(allEdges, 0);

    forAll(faceEdges_, faceI)
    {
        const labelList& faceEdges = faceEdges_[faceI];

        if (faceEdges.empty())
        {
            continue;
        }

        // Check consistency of face-edge-points as well
        edgeList eList = faces_[faceI].edges();

        forAll(faceEdges,edgeI)
        {
            nEdgeFaces[faceEdges[edgeI]]++;

            // Check if this edge actually belongs to this face
            bool found = false;
            const edge& edgeToCheck = edges_[faceEdges[edgeI]];

            forAll(eList, edgeII)
            {
                if (edgeToCheck == eList[edgeII])
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                Pout << nl << nl << "Edge: " << faceEdges[edgeI]
                     << ": " << edgeToCheck << nl
                     << "was not found in face: " << faceI
                     << ": " << faces_[faceI] << nl
                     << "faceEdges: " << faceEdges
                     << endl;

                nFailedChecks++;

                ConnectivityWarning()
                     << nl << "Edge-Face connectivity is inconsistent."
                     << endl;
            }
        }
    }

    label nInternalEdges = 0;
    labelList patchInfo(boundaryMesh().size(), 0);

    forAll(edgeFaces_, edgeI)
    {
        const labelList& edgeFaces = edgeFaces_[edgeI];

        if (edgeFaces.empty())
        {
            continue;
        }

        if (edgeFaces.size() != nEdgeFaces[edgeI])
        {
            Pout << nl << nl << "Edge: " << edgeI
                 << ": edgeFaces: " << edgeFaces << endl;

            nFailedChecks++;

            ConnectivityWarning()
                << nl << "Edge-Face connectivity is inconsistent."
                << endl;
        }

        label nBF = 0;

        // Check if this edge belongs to faceEdges for each face
        forAll(edgeFaces, faceI)
        {
            if (findIndex(faceEdges_[edgeFaces[faceI]], edgeI) == -1)
            {
                Pout << nl << nl << "Edge: " << edgeI << ": " << edges_[edgeI]
                     << ", edgeFaces: " << edgeFaces << nl
                     << "was not found in faceEdges of face: "
                     << edgeFaces[faceI] << ": " << faces_[edgeFaces[faceI]]
                     << nl << "faceEdges: " << faceEdges_[edgeFaces[faceI]]
                     << endl;

                nFailedChecks++;

                ConnectivityWarning()
                    << nl << "Edge-Face connectivity is inconsistent."
                    << endl;
            }

            if (neighbour_[edgeFaces[faceI]] == -1)
            {
                nBF++;
            }
        }

        if (nBF == 0)
        {
            nInternalEdges++;

            // Check if this edge is actually internal.
            if (whichEdgePatch(edgeI) >= 0)
            {
                Pout << "Edge: " << edgeI
                     << ": " << edges_[edgeI] << " is internal, "
                     << " but patch is specified as: "
                     << whichEdgePatch(edgeI)
                     << endl;

                nFailedChecks++;
            }
        }
        else
        {
            label patchID = whichEdgePatch(edgeI);

            // Check if this edge is actually on a boundary.
            if (patchID < 0)
            {
                Pout << "Edge: " << edgeI
                     << ": " << edges_[edgeI]
                     << " is on a boundary, but patch is specified as: "
                     << patchID << endl;

                nFailedChecks++;
            }
            else
            {
                patchInfo[patchID]++;
            }

            if (nBF > 2)
            {
                Pout << "Edge: " << edgeI
                     << ": " << edges_[edgeI]
                     << " has " << nBF
                     << " boundary faces connected to it." << nl
                     << " Pinched manifolds are not allowed."
                     << endl;

                nFailedChecks++;
            }
        }
    }

    if (nInternalEdges != nInternalEdges_)
    {
        Pout << nl << "Internal edge-count is inconsistent." << nl
             << " Counted internal edges: " << nInternalEdges
             << " Actual count: " << nInternalEdges_ << endl;

        nFailedChecks++;
    }

    forAll(patchInfo, patchI)
    {
        if (patchInfo[patchI] != edgePatchSizes_[patchI])
        {
            Pout << "Patch-count is inconsistent." << nl
                 << " Patch: " << patchI
                 << " Counted edges: " << patchInfo[patchI]
                 << " Actual count: " << edgePatchSizes_[patchI] << endl;

            nFailedChecks++;
        }
    }

    // Check added edge patches to ensure that it is consistent
    forAllConstIter(Map<label>, addedEdgePatches_, aepIter)
    {
        label key = aepIter.key();
        label patch = aepIter();

        label nBF = 0;
        const labelList& edgeFaces = edgeFaces_[key];

        // Check if any faces on boundaries
        forAll(edgeFaces, faceI)
        {
            if (neighbour_[edgeFaces[faceI]] == -1)
            {
                nBF++;
            }
        }

        if ((patch < 0) && (nBF > 0))
        {
            Pout << nl << nl << "Edge: " << key
                 << ", edgeFaces: " << edgeFaces
                 << " is internal, but contains boundary faces."
                 << endl;

            nFailedChecks++;
        }

        if ((patch >= 0) && (nBF != 2))
        {
            Pout << nl << nl << "Edge: " << key
                 << ", edgeFaces: " << edgeFaces
                 << " is on a boundary patch, but doesn't contain"
                 << " two boundary faces."
                 << endl;

            nFailedChecks++;
        }
    }

    Info << "Done." << endl;

    if (!twoDMesh_)
    {
        Info << "Checking point-edge connectivity...";

        label allPoints = points_.size();
        List<labelHashSet> hlPointEdges(allPoints);

        forAll(edges_, edgeI)
        {
            if (edgeFaces_[edgeI].size())
            {
                hlPointEdges[edges_[edgeI][0]].insert(edgeI);
                hlPointEdges[edges_[edgeI][1]].insert(edgeI);
            }
        }

        forAll(pointEdges_, pointI)
        {
            const labelList& pointEdges = pointEdges_[pointI];

            if (pointEdges.empty())
            {
                continue;
            }

            forAll(pointEdges, edgeI)
            {
                if (!hlPointEdges[pointI].found(pointEdges[edgeI]))
                {
                    Pout << nl << nl << "Point: " << pointI << nl
                         << "pointEdges: " << pointEdges << nl
                         << "hlPointEdges: " << hlPointEdges[pointI]
                         << endl;

                    nFailedChecks++;

                    ConnectivityWarning()
                        << nl << "Point-Edge connectivity is inconsistent."
                        << endl;
                }
            }

            // Do a size check as well
            if
            (
                hlPointEdges[pointI].size() != pointEdges.size() ||
                pointEdges.size() == 1
            )
            {
                Pout << nl << nl << "Point: " << pointI << nl
                     << "pointEdges: " << pointEdges << nl
                     << "hlPointEdges: " << hlPointEdges[pointI]
                     << endl;

                nFailedChecks++;

                ConnectivityWarning()
                    << "Size inconsistency."
                    << nl << "Point-Edge connectivity is inconsistent."
                    << endl;
            }
        }

        Info << "Done." << endl;

        Info << "Checking edge-points connectivity...";

        label otherPoint = -1, nextPoint = -1;

        forAll(edgePoints_, edgeI)
        {
            // Do a preliminary size check
            const labelList& edgePoints = edgePoints_[edgeI];
            const labelList& edgeFaces = edgeFaces_[edgeI];

            if (edgeFaces.empty())
            {
                continue;
            }

            if (edgePoints.size() != edgeFaces.size())
            {
                Pout << nl << nl
                     << "Edge: " << edgeI
                     << " " << edges_[edgeI] << endl;

                Pout << "edgeFaces: " << edgeFaces << endl;
                forAll(edgeFaces, faceI)
                {
                    Info << edgeFaces[faceI] << ": "
                         << faces_[edgeFaces[faceI]]
                         << endl;
                }

                Pout << "edgePoints: " << edgePoints << endl;

                nFailedChecks++;

                ConnectivityWarning()
                    << nl << "Edge-Points connectivity is inconsistent."
                    << endl;
            }

            // Now check to see that both lists are consistent.
            const edge& edgeToCheck = edges_[edgeI];

            forAll(edgeFaces, faceI)
            {
                const face& faceToCheck = faces_[edgeFaces[faceI]];

                meshOps::findIsolatedPoint
                (
                    faceToCheck,
                    edgeToCheck,
                    otherPoint,
                    nextPoint
                );

                if (findIndex(edgePoints, otherPoint) == -1)
                {
                    Pout << nl << nl
                         << "Edge: " << edgeI
                         << " " << edges_[edgeI] << endl;

                    Pout << "edgeFaces: " << edgeFaces << endl;
                    forAll(edgeFaces, faceI)
                    {
                        Info << edgeFaces[faceI] << ": "
                             << faces_[edgeFaces[faceI]]
                             << endl;
                    }

                    Pout << "edgePoints: " << edgePoints << endl;

                    nFailedChecks++;

                    ConnectivityWarning()
                        << nl << "Edge-Points connectivity is inconsistent."
                        << endl;
                }
            }
        }

        Info << "Done." << endl;
    }

    Info << "Checking cell-point connectivity...";

    // Loop through all cells and construct cell-to-node
    label cIndex = 0;
    label allCells = cells_.size();
    labelList cellIndex(allCells);
    List<labelHashSet> cellToNode(allCells);

    forAll(cells_, cellI)
    {
        const cell& thisCell = cells_[cellI];

        if (thisCell.empty())
        {
            continue;
        }

        cellIndex[cIndex] = cellI;

        forAll(thisCell, faceI)
        {
            const labelList& fEdges = faceEdges_[thisCell[faceI]];

            forAll(fEdges, edgeI)
            {
                const edge& thisEdge = edges_[fEdges[edgeI]];

                if (!cellToNode[cIndex].found(thisEdge[0]))
                {
                    cellToNode[cIndex].insert(thisEdge[0]);
                }

                if (!cellToNode[cIndex].found(thisEdge[1]))
                {
                    cellToNode[cIndex].insert(thisEdge[1]);
                }
            }
        }

        cIndex++;
    }

    // Resize the lists
    cellIndex.setSize(cIndex);
    cellToNode.setSize(cIndex);

    // Preliminary check for size
    forAll(cellToNode, cellI)
    {
        if
        (
            (cellToNode[cellI].size() != 6 && twoDMesh_) ||
            (cellToNode[cellI].size() != 4 && !twoDMesh_)
        )
        {
            Pout << nl << "Warning: Cell: "
                 << cellIndex[cellI] << " is inconsistent. "
                 << endl;

            const cell& failedCell = cells_[cellIndex[cellI]];

            Info << "Cell faces: " << failedCell << endl;

            forAll(failedCell, faceI)
            {
                Info << "\tFace: " << failedCell[faceI]
                     << " :: " << faces_[failedCell[faceI]]
                     << endl;

                const labelList& fEdges = faceEdges_[failedCell[faceI]];

                forAll(fEdges, edgeI)
                {
                    Info << "\t\tEdge: " << fEdges[edgeI]
                         << " :: " << edges_[fEdges[edgeI]]
                         << endl;
                }
            }

            nFailedChecks++;
        }
    }

    Info << "Done." << endl;

    reduce(nFailedChecks, orOp<bool>());

    if (nFailedChecks)
    {
        FatalErrorIn
        (
            "void dynamicTopoFvMesh::checkConnectivity"
            "(const label maxErrors) const"
        )
            << nFailedChecks << " failures were found in connectivity."
            << abort(FatalError);
    }
}


// Utility method to check the quality
// of a triangular face after bisection.
//  - Returns 'true' if the bisection in NOT feasible.
bool dynamicTopoFvMesh::checkBisection
(
    const label fIndex,
    const label bFaceIndex,
    bool forceOp
) const
{
    scalar bisectionQuality = GREAT, minArea = GREAT;

    label commonEdge = -1;

    // Find the common edge index
    meshOps::findCommonEdge
    (
        bFaceIndex,
        fIndex,
        faceEdges_,
        commonEdge
    );

    // Fetch the edge
    const edge& checkEdge = edges_[commonEdge];
    const face& checkFace = faces_[bFaceIndex];

    // Compute old / new mid-points
    point mpOld =
    (
        linePointRef
        (
            oldPoints_[checkEdge.start()],
            oldPoints_[checkEdge.end()]
        ).centre()
    );

    point mpNew =
    (
        linePointRef
        (
            points_[checkEdge.start()],
            points_[checkEdge.end()]
        ).centre()
    );

    // Find the isolated point on the face
    label iPoint = -1, nPoint = -1;

    meshOps::findIsolatedPoint
    (
        checkFace,
        checkEdge,
        iPoint,
        nPoint
    );

    // Find the other point
    label oPoint =
    (
        (nPoint == checkEdge.start()) ?
        checkEdge.end() : checkEdge.start()
    );

    // Configure old / new triangle faces
    FixedList<FixedList<point, 3>, 2> tfNew, tfOld;

    tfNew[0][0] = mpNew;
    tfNew[0][1] = points_[iPoint];
    tfNew[0][2] = points_[nPoint];

    tfOld[0][0] = mpOld;
    tfOld[0][1] = oldPoints_[iPoint];
    tfOld[0][2] = oldPoints_[nPoint];

    tfNew[1][0] = points_[oPoint];
    tfNew[1][1] = points_[iPoint];
    tfNew[1][2] = mpNew;

    tfOld[1][0] = oldPoints_[oPoint];
    tfOld[1][1] = oldPoints_[iPoint];
    tfOld[1][2] = mpOld;

    // Assume XY plane here
    vector n = vector(0,0,1);

    forAll(tfNew, fI)
    {
        // Configure triangles
        triPointRef tprNew(tfNew[fI][0], tfNew[fI][1], tfNew[fI][2]);
        triPointRef tprOld(tfOld[fI][0], tfOld[fI][1], tfOld[fI][2]);

        scalar tQuality =
        (
            tprNew.quality() *
            (
                Foam::sign
                (
                    tprNew.normal() &
                    ((tprNew.centre() & n) * n)
                )
            )
        );

        scalar oldArea =
        (
            mag(tprOld.normal()) *
            (
                Foam::sign
                (
                    tprOld.normal() &
                    ((tprOld.centre() & n) * n)
                )
            )
        );

        // Update statistics
        minArea = Foam::min(minArea, oldArea);
        bisectionQuality = Foam::min(bisectionQuality, tQuality);
    }

    // Final quality check
    if (bisectionQuality < sliverThreshold_ && !forceOp)
    {
        return true;
    }

    // Negative quality is a no-no
    if (bisectionQuality < 0.0)
    {
        return true;
    }

    // Negative old-area is also a no-no
    if (minArea < 0.0)
    {
        return true;
    }

    // No problems, so a bisection is feasible.
    return false;
}


// Utility method to check whether the cell given by 'cellIndex' will yield
// a valid cell when 'pointIndex' is moved to 'newPoint'.
//  - The routine performs metric-based checks.
//  - Returns 'true' if the collapse in NOT feasible, and
//    makes entries in cellsChecked to avoid repetitive checks.
bool dynamicTopoFvMesh::checkCollapse
(
    const labelList& triFaces,
    const FixedList<label,2>& c0BdyIndex,
    const FixedList<label,2>& c1BdyIndex,
    const FixedList<label,2>& pointIndex,
    const FixedList<point,2>& newPoint,
    const FixedList<point,2>& oldPoint,
    scalar& collapseQuality,
    const bool checkNeighbour,
    bool forceOp
) const
{
    // Reset input
    collapseQuality = GREAT;
    scalar minArea = GREAT;

    forAll(triFaces, indexI)
    {
        if
        (
            (triFaces[indexI] == c0BdyIndex[0])
         || (triFaces[indexI] == c0BdyIndex[1])
        )
        {
            continue;
        }

        if (checkNeighbour)
        {
            if
            (
                (triFaces[indexI] == c1BdyIndex[0])
             || (triFaces[indexI] == c1BdyIndex[1])
            )
            {
                continue;
            }
        }

        const face& checkFace = faces_[triFaces[indexI]];

        // Configure a triangle face
        FixedList<point, 3> tFNew(vector::zero);
        FixedList<point, 3> tFOld(vector::zero);

        // Make necessary replacements
        forAll(checkFace, pointI)
        {
            tFNew[pointI] = points_[checkFace[pointI]];
            tFOld[pointI] = oldPoints_[checkFace[pointI]];

            if (checkFace[pointI] == pointIndex[0])
            {
                tFNew[pointI] = newPoint[0];
                tFOld[pointI] = oldPoint[0];
            }

            if (checkFace[pointI] == pointIndex[1])
            {
                tFNew[pointI] = newPoint[1];
                tFOld[pointI] = oldPoint[1];
            }
        }

        // Configure triangles
        triPointRef tprNew(tFNew[0], tFNew[1], tFNew[2]);
        triPointRef tprOld(tFOld[0], tFOld[1], tFOld[2]);

        // Assume XY plane here
        vector n = vector(0,0,1);

        // Compute the quality.
        // Assume centre-plane passes through origin
        scalar tQuality =
        (
            tprNew.quality() *
            (
                Foam::sign
                (
                    tprNew.normal() &
                    ((tprNew.centre() & n) * n)
                )
            )
        );

        scalar oldArea =
        (
            mag(tprOld.normal()) *
            (
                Foam::sign
                (
                    tprOld.normal() &
                    ((tprOld.centre() & n) * n)
                )
            )
        );

        // Update statistics
        minArea = Foam::min(minArea, oldArea);
        collapseQuality = Foam::min(collapseQuality, tQuality);
    }

    // Final quality check
    if (collapseQuality < sliverThreshold_ && !forceOp)
    {
        return true;
    }

    // Negative quality is a no-no
    if (collapseQuality < 0.0)
    {
        return true;
    }

    // Negative old-area is also a no-no
    if (minArea < 0.0)
    {
        return true;
    }

    // No problems, so a collapse is feasible.
    return false;
}


// Utility method to check whether the cell given by 'cellIndex' will yield
// a valid cell when 'pointIndex' is moved to 'newPoint'.
//  - The routine performs metric-based checks.
//  - Returns 'true' if the collapse in NOT feasible, and
//    makes entries in cellsChecked to avoid repetitive checks.
bool dynamicTopoFvMesh::checkCollapse
(
    const point& newPoint,
    const point& oldPoint,
    const label pointIndex,
    const label cellIndex,
    labelHashSet& cellsChecked,
    scalar& collapseQuality,
    bool forceOp
) const
{
    label faceIndex = -1;
    scalar cQuality = 0.0, oldVolume = 0.0;
    const cell& cellToCheck = cells_[cellIndex];

    // Look for a face that doesn't contain 'pointIndex'
    forAll(cellToCheck, faceI)
    {
        const face& currFace = faces_[cellToCheck[faceI]];

        if (currFace.which(pointIndex) < 0)
        {
            faceIndex = cellToCheck[faceI];
            break;
        }
    }

    // Compute cell-volume
    const face& faceToCheck = faces_[faceIndex];

    if (owner_[faceIndex] == cellIndex)
    {
        cQuality =
        (
            tetMetric_
            (
                points_[faceToCheck[2]],
                points_[faceToCheck[1]],
                points_[faceToCheck[0]],
                newPoint
            )
        );

        oldVolume =
        (
            tetPointRef
            (
                oldPoints_[faceToCheck[2]],
                oldPoints_[faceToCheck[1]],
                oldPoints_[faceToCheck[0]],
                oldPoint
            ).mag()
        );
    }
    else
    {
        cQuality =
        (
            tetMetric_
            (
                points_[faceToCheck[0]],
                points_[faceToCheck[1]],
                points_[faceToCheck[2]],
                newPoint
            )
        );

        oldVolume =
        (
            tetPointRef
            (
                oldPoints_[faceToCheck[0]],
                oldPoints_[faceToCheck[1]],
                oldPoints_[faceToCheck[2]],
                oldPoint
            ).mag()
        );
    }

    // Final quality check
    if (cQuality < sliverThreshold_ && !forceOp)
    {
        if (debug > 3)
        {
            InfoIn
            (
                "\n\n"
                "bool dynamicTopoFvMesh::checkCollapse\n"
                "(\n"
                "    const point& newPoint,\n"
                "    const point& oldPoint,\n"
                "    const label pointIndex,\n"
                "    const label cellIndex,\n"
                "    labelHashSet& cellsChecked,\n"
                "    scalar& collapseQuality,\n"
                "    bool forceOp\n"
                ") const\n"
            )
                << "\nCollapsing cell: " << cellIndex
                << " containing points:\n"
                << faceToCheck[0] << "," << faceToCheck[1] << ","
                << faceToCheck[2] << "," << pointIndex << nl
                << "will yield a quality of: " << cQuality
                << ", when " << pointIndex
                << " is moved to location: " << nl
                << newPoint
                << endl;
        }

        return true;
    }

    // Negative quality is a no-no
    if (cQuality < 0.0)
    {
        if (forceOp)
        {
            InfoIn
            (
                "\n\n"
                "bool dynamicTopoFvMesh::checkCollapse\n"
                "(\n"
                "    const point& newPoint,\n"
                "    const point& oldPoint,\n"
                "    const label pointIndex,\n"
                "    const label cellIndex,\n"
                "    labelHashSet& cellsChecked,\n"
                "    scalar& collapseQuality,\n"
                "    bool forceOp\n"
                ") const\n"
            )
                << "\nCollapsing cell: " << cellIndex
                << " containing points:\n"
                << faceToCheck[0] << "," << faceToCheck[1] << ","
                << faceToCheck[2] << "," << pointIndex << nl
                << "will yield a quality of: " << cQuality
                << ", when " << pointIndex
                << " is moved to location: " << nl
                << newPoint << nl
                << "Operation cannot be forced."
                << endl;
        }

        return true;
    }

    // Negative old-volume is also a no-no
    if (oldVolume < 0.0)
    {
        if (forceOp)
        {
            InfoIn
            (
                "\n\n"
                "bool dynamicTopoFvMesh::checkCollapse\n"
                "(\n"
                "    const point& newPoint,\n"
                "    const point& oldPoint,\n"
                "    const label pointIndex,\n"
                "    const label cellIndex,\n"
                "    labelHashSet& cellsChecked,\n"
                "    scalar& collapseQuality,\n"
                "    bool forceOp\n"
                ") const\n"
            )
                << "\nCollapsing cell: " << cellIndex
                << " containing points:\n"
                << faceToCheck[0] << "," << faceToCheck[1] << ","
                << faceToCheck[2] << "," << pointIndex << nl
                << "will yield an old-volume of: " << oldVolume
                << ", when " << pointIndex
                << " is moved to location: " << nl
                << oldPoint << nl
                << "Operation cannot be forced."
                << endl;
        }

        return true;
    }

    // No problems, so a collapse is feasible
    cellsChecked.insert(cellIndex);

    // Update input quality
    collapseQuality = Foam::min(collapseQuality, cQuality);

    return false;
}


} // End namespace Foam

// ************************************************************************* //
