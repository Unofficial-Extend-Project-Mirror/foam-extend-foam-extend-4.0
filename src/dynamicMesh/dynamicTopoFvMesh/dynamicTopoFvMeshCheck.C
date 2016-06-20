/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    dynamicTopoFvMesh

Description
    Functions specific to connectivity checking and debugging

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "dynamicTopoFvMesh.H"

#include "IOmanip.H"
#include "volFields.H"
#include "triPointRef.H"
#include "tetPointRef.H"
#include "coupledInfo.H"

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

        // Skip hexahedral cells
        if (cellToCheck.size() == 6)
        {
            cQuality = 1.0;
            meanQuality += cQuality;

            // Update min / max
            maxQuality = Foam::max(cQuality, maxQuality);
            minQuality = Foam::min(cQuality, minQuality);

            nCells++;

            continue;
        }

        if (is2D())
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
        if (minQuality < 0.0)
        {
            Pout<< nl
                << " Warning: Minimum cell quality is: " << minQuality
                << " at cell: " << minCell
                << endl;
        }

        // Reduce statistics across processors.
        reduce(minQuality, minOp<scalar>());
        reduce(maxQuality, maxOp<scalar>());
        reduce(meanQuality, sumOp<scalar>());
        reduce(nCells, sumOp<label>());

        Info<< nl
            << " ~~~ Mesh Quality Statistics ~~~ " << nl
            << " Min: " << minQuality << nl
            << " Max: " << maxQuality << nl
            << " Mean: " << meanQuality/nCells << nl
            << " Cells: " << nCells << nl
            << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << nl
            << endl;
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
//  - If nProcCurves is provided, the variable is incremented
//    if the edge is processor-coupled
bool dynamicTopoFvMesh::checkBoundingCurve
(
    const label eIndex,
    const bool overRidePurityCheck,
    label* nProcCurves
) const
{
    // Internal edges don't count
    label edgePatch = -1;

    // If this entity was deleted, skip it.
    if (edgeFaces_[eIndex].empty())
    {
        // Return true so that swap3DEdges skips this edge.
        return true;
    }

    // Check if two boundary faces lie on different face-patches
    bool procCoupled = false;
    FixedList<label, 2> fPatches(-1);
    FixedList<vector, 2> fNorm(vector::zero);

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

        // Explicit check for processor edges (both 2D and 3D)
        if (processorCoupledEntity(eIndex, false, true))
        {
            // Increment nProcCurves
            if (nProcCurves)
            {
                (*nProcCurves)++;
            }

            // Check for pure processor edge, and if not,
            // fetch boundary patch labels / normals
            if
            (
                processorCoupledEntity
                (
                    eIndex,
                    false,
                    true,
                    true,
                    &fPatches,
                    &fNorm
                )
            )
            {
                // 'Pure' processor coupled edges don't count
                return false;
            }
            else
            if (!overRidePurityCheck)
            {
                // This edge lies between a processor and physical patch,
                //  - This a bounding curve (unless an override is requested)
                //  - An override is warranted for 2-2 swaps on impure edges,
                //    which is typically requested by swap3DEdges.
                return true;
            }

            // Specify that the edge is procCoupled
            procCoupled = true;
        }
    }

    if (procCoupled)
    {
        // Normalize patch normals from coupled check
        fNorm[0] /= mag(fNorm[0]) + VSMALL;
        fNorm[1] /= mag(fNorm[1]) + VSMALL;
    }
    else
    {
        // Fetch patch indices / normals
        const labelList& eFaces = edgeFaces_[eIndex];

        label fPatch = -1, count = 0;

        forAll(eFaces, faceI)
        {
            if ((fPatch = whichPatch(eFaces[faceI])) > -1)
            {
                // Obtain the normal.
                fNorm[count] = faces_[eFaces[faceI]].normal(points_);

                // Normalize it.
                fNorm[count] /= mag(fNorm[count]) + VSMALL;

                // Note patch index
                fPatches[count] = fPatch;

                count++;

                if (count == 2)
                {
                    break;
                }
            }
        }
    }

    // Check for legitimate patches
    if (fPatches[0] < 0 || fPatches[1] < 0)
    {
        const labelList& eFaces = edgeFaces_[eIndex];

        forAll(eFaces, faceI)
        {
            Pout<< " Face: " << eFaces[faceI]
                << " :: " << faces_[eFaces[faceI]]
                << " Patch: " << whichPatch(eFaces[faceI]) << nl;
        }

        label epI = whichEdgePatch(eIndex);

        FatalErrorIn
        (
            "bool dynamicTopoFvMesh::checkBoundingCurve"
            "(const label, const bool) const"
        )
            << " Edge: " << eIndex << ":: " << edges_[eIndex]
            << " Patch: "
            << (epI < 0 ? "Internal" : boundaryMesh()[epI].name())
            << " edgeFaces: " << eFaces << nl
            << " expected 2 boundary patches." << nl
            << " fPatches[0]: " << fPatches[0] << nl
            << " fPatches[1]: " << fPatches[1] << nl
            << " fNorm[0]: " << fNorm[0] << nl
            << " fNorm[1]: " << fNorm[1] << nl
            << " coupledModification: " << coupledModification_
            << abort(FatalError);
    }

    scalar deviation = (fNorm[0] & fNorm[1]);

    // Check if the swap-curvature is too high
    if (mag(deviation) < swapDeviation_)
    {
        return true;
    }

    // Check if the edge borders two different patches
    if (fPatches[0] != fPatches[1])
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
            Pout<< nl << nl
                << " eIndex: " << eIndex
                << " minQuality: " << minQuality
                << " newQuality: " << Q[checkIndex][0][m[checkIndex]-1]
                << endl;
        }
    }

    if (coupledModification_)
    {
        // Only locally coupled indices require checks
        if (locallyCoupledEntity(eIndex))
        {
            // Check the quality of the slave edge as well.
            label sIndex = -1;

            // Loop through masterToSlave and determine the slave index.
            forAll(patchCoupling_, patchI)
            {
                if (patchCoupling_(patchI))
                {
                    const label edgeEnum  = coupleMap::EDGE;
                    const coupleMap& cMap = patchCoupling_[patchI].map();

                    if ((sIndex = cMap.findSlave(edgeEnum, eIndex)) > -1)
                    {
                        break;
                    }
                }
            }

            if (sIndex == -1)
            {
                FatalErrorIn
                (
                    "bool dynamicTopoFvMesh::checkQuality\n"
                    "(\n"
                    "    const label eIndex,\n"
                    "    const labelList& m,\n"
                    "    const PtrList<scalarListList>& Q,\n"
                    "    const scalar minQuality,\n"
                    "    const label checkIndex\n"
                    ") const\n"
                )
                    << "Coupled maps were improperly specified." << nl
                    << " Slave index not found for: " << nl
                    << " Edge: " << eIndex << nl
                    << abort(FatalError);
            }

            // Turn off switch temporarily.
            unsetCoupledModification();

            // Recursively call for the slave edge.
            myResult =
            (
                myResult && checkQuality(sIndex, m, Q, minQuality, 1)
            );

            // Turn it back on.
            setCoupledModification();
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
    Pout<< "m: " << m[checkIndex] << endl;

    // Print out Q
    Pout<< "===" << nl
        << " Q " << nl
        << "===" << endl;

    Pout<< "   ";

    for (label j = 0; j < m[checkIndex]; j++)
    {
        Pout<< setw(12) << j;
    }

    Pout<< nl;

    for (label j = 0; j < m[checkIndex]; j++)
    {
        Pout<< "-------------";
    }

    Pout<< nl;

    for (label i = 0; i < (m[checkIndex]-2); i++)
    {
        Pout<< i << ": ";

        for (label j = 0; j < m[checkIndex]; j++)
        {
            Pout<< setw(12) << Q[checkIndex][i][j];
        }

        Pout<< nl;
    }

    // Print out K
    Pout<< "===" << nl
        << " K " << nl
        << "===" << endl;

    Pout<< "   ";

    for (label j = 0; j < m[checkIndex]; j++)
    {
        Pout<< setw(12) << j;
    }

    Pout<< nl;

    for (label i = 0; i < (m[checkIndex]-2); i++)
    {
        Pout<< i << ": ";

        for (label j = 0; j < m[checkIndex]; j++)
        {
            Pout<< setw(12) << K[checkIndex][i][j];
        }

        Pout<< nl;
    }

    Pout<< endl;
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
    scalar oldTetVol = 0.0, newTetVol = 0.0;

    const edge& edgeToCheck = edges_[eIndex];

    for (label i = 0; i < (m-2); i++)
    {
        // Compute volume for the upper-half
        newTetVol =
        (
            tetPointRef
            (
                points_[hullVertices[triangulations[0][i]]],
                points_[hullVertices[triangulations[1][i]]],
                points_[hullVertices[triangulations[2][i]]],
                points_[edgeToCheck[0]]
            ).mag()
        );

        oldTetVol =
        (
            tetPointRef
            (
                oldPoints_[hullVertices[triangulations[0][i]]],
                oldPoints_[hullVertices[triangulations[1][i]]],
                oldPoints_[hullVertices[triangulations[2][i]]],
                oldPoints_[edgeToCheck[0]]
            ).mag()
        );

        if (oldTetVol < 0.0 || (mag(oldTetVol) < mag(0.1 * newTetVol)))
        {
            if (debug > 2)
            {
                Pout<< " Swap sequence leads to bad old-volumes." << nl
                    << " Edge: " << edgeToCheck << nl
                    << " using Points: " << nl
                    << oldPoints_[hullVertices[triangulations[0][i]]] << nl
                    << oldPoints_[hullVertices[triangulations[1][i]]] << nl
                    << oldPoints_[hullVertices[triangulations[2][i]]] << nl
                    << oldPoints_[edgeToCheck[0]] << nl
                    << " Old Volume: " << oldTetVol << nl
                    << " New Volume: " << newTetVol << nl
                    << endl;
            }

            return true;
        }

        newTetVol =
        (
            tetPointRef
            (
                points_[hullVertices[triangulations[2][i]]],
                points_[hullVertices[triangulations[1][i]]],
                points_[hullVertices[triangulations[0][i]]],
                points_[edgeToCheck[1]]
            ).mag()
        );

        oldTetVol =
        (
            tetPointRef
            (
                oldPoints_[hullVertices[triangulations[2][i]]],
                oldPoints_[hullVertices[triangulations[1][i]]],
                oldPoints_[hullVertices[triangulations[0][i]]],
                oldPoints_[edgeToCheck[1]]
            ).mag()
        );

        if (oldTetVol < 0.0 || (mag(oldTetVol) < mag(0.1 * newTetVol)))
        {
            if (debug > 2)
            {
                Pout<< " Swap sequence leads to bad old-volumes." << nl
                    << " Edge: " << edgeToCheck << nl
                    << " using Points: " << nl
                    << oldPoints_[hullVertices[triangulations[2][i]]] << nl
                    << oldPoints_[hullVertices[triangulations[1][i]]] << nl
                    << oldPoints_[hullVertices[triangulations[0][i]]] << nl
                    << oldPoints_[edgeToCheck[1]] << nl
                    << " Old Volume: " << oldTetVol << nl
                    << " New Volume: " << newTetVol << nl
                    << endl;
            }

            return true;
        }
    }

    return false;
}


// Write out connectivity for an edge
void dynamicTopoFvMesh::writeEdgeConnectivity
(
    const label eIndex
) const
{
    // Write out edge
    writeVTK("Edge_" + Foam::name(eIndex), eIndex, 1, false, true);

    const labelList& eFaces = edgeFaces_[eIndex];

    // Write out edge faces
    writeVTK
    (
        "EdgeFaces_"
      + Foam::name(eIndex)
      + '_'
      + Foam::name(Pstream::myProcNo()),
        eFaces,
        2, false, true
    );

    dynamicLabelList edgeCells(10);

    forAll(eFaces, faceI)
    {
        label pIdx = whichPatch(eFaces[faceI]);

        word pName((pIdx < 0) ? "Internal" : boundaryMesh()[pIdx].name());

        Pout<< " Face: " << eFaces[faceI]
            << " :: " << faces_[eFaces[faceI]]
            << " Patch: " << pName
            << " Proc: " << Pstream::myProcNo() << nl;

        label own = owner_[eFaces[faceI]];
        label nei = neighbour_[eFaces[faceI]];

        if (findIndex(edgeCells, own) == -1)
        {
            edgeCells.append(own);
        }

        if (nei == -1)
        {
            continue;
        }

        if (findIndex(edgeCells, nei) == -1)
        {
            edgeCells.append(nei);
        }
    }

    // Write out cells connected to edge
    writeVTK
    (
        "EdgeCells_"
      + Foam::name(eIndex)
      + '_'
      + Foam::name(Pstream::myProcNo()),
        edgeCells,
        3, false, true
    );

    if (is2D())
    {
        return;
    }

    // Check processors for coupling
    forAll(procIndices_, pI)
    {
        // Fetch reference to subMesh
        const coupleMap& cMap = recvMeshes_[pI].map();
        const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

        label sI = -1;

        if ((sI = cMap.findSlave(coupleMap::EDGE, eIndex)) == -1)
        {
            continue;
        }

        const edge& slaveEdge = mesh.edges_[sI];
        const labelList& seFaces = mesh.edgeFaces_[sI];

        edge cE
        (
            cMap.findMaster(coupleMap::POINT, slaveEdge[0]),
            cMap.findMaster(coupleMap::POINT, slaveEdge[1])
        );

        Pout<< " >> Edge: " << sI << "::" << slaveEdge
            << " mapped: " << cE << nl;

        mesh.writeVTK
        (
            "EdgeFaces_"
          + Foam::name(eIndex)
          + '_'
          + Foam::name(procIndices_[pI]),
            seFaces,
            2, false, true
        );

        // Clear existing list
        edgeCells.clear();

        forAll(seFaces, faceI)
        {
            label pIdx = mesh.whichPatch(seFaces[faceI]);

            word pName
            (
                (pIdx < 0) ?
                "Internal" :
                mesh.boundaryMesh()[pIdx].name()
            );

            Pout<< " Face: " << seFaces[faceI]
                << " :: " << mesh.faces_[seFaces[faceI]]
                << " Patch: " << pName
                << " Proc: " << procIndices_[pI] << nl;

            label own = mesh.owner_[seFaces[faceI]];
            label nei = mesh.neighbour_[seFaces[faceI]];

            if (findIndex(edgeCells, own) == -1)
            {
                edgeCells.append(own);
            }

            if (nei == -1)
            {
                continue;
            }

            if (findIndex(edgeCells, nei) == -1)
            {
                edgeCells.append(nei);
            }
        }

        mesh.writeVTK
        (
            "EdgeCells_"
          + Foam::name(eIndex)
          + '_'
          + Foam::name(procIndices_[pI]),
            edgeCells,
            3, false, true
        );
    }
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
    const bool useOldPoints,
    const UList<scalar>& scalField,
    const UList<label>& lablField
) const
{
    // Check if spatial bounding box has been specified
    const dictionary& meshSubDict = dict_.subDict("dynamicTopoFvMesh");

    labelList entityList;

    if (meshSubDict.found("spatialDebug") && !useOldConnectivity)
    {
        // Read the bounding box
        boundBox bb
        (
            meshSubDict.subDict("spatialDebug").lookup("debugBoundBox")
        );

        dynamicLabelList cSubList(10);

        forAll(cList, cellI)
        {
            label index = cList[cellI];

            if (index < 0)
            {
                continue;
            }

            point containPoint(vector::zero);

            switch (primitiveType)
            {
                // Are we looking at points?
                case 0:
                {
                    containPoint = points_[index];
                    break;
                }

                // Are we looking at edges?
                case 1:
                {
                    containPoint = edges_[index].centre(points_);
                    break;
                }

                // Are we looking at faces?
                case 2:
                {
                    containPoint = faces_[index].centre(points_);
                    break;
                }

                // Are we looking at cells?
                case 3:
                {
                    scalar volume = 0.0;

                    // Compute centre
                    meshOps::cellCentreAndVolume
                    (
                        index,
                        points_,
                        faces_,
                        cells_,
                        owner_,
                        containPoint,
                        volume
                    );

                    break;
                }
            }

            // Is the point of interest?
            if (bb.contains(containPoint))
            {
                cSubList.append(index);
            }
        }

        // If nothing is present, don't write out anything
        if (cSubList.empty())
        {
            return;
        }
        else
        {
            // Take over contents
            entityList = cSubList;
        }
    }
    else
    {
        // Conventional output
        entityList = cList;
    }

    if (useOldPoints)
    {
        if (useOldConnectivity)
        {
            // Use points from polyMesh
            meshOps::writeVTK
            (
                (*this),
                name,
                entityList,
                primitiveType,
                polyMesh::points(),
                polyMesh::edges(),
                polyMesh::faces(),
                polyMesh::cells(),
                polyMesh::faceOwner(),
                scalField,
                lablField
            );
        }
        else
        {
            meshOps::writeVTK
            (
                (*this),
                name,
                entityList,
                primitiveType,
                oldPoints_,
                edges_,
                faces_,
                cells_,
                owner_,
                scalField,
                lablField
            );
        }
    }
    else
    {
        meshOps::writeVTK
        (
            (*this),
            name,
            entityList,
            primitiveType,
            points_,
            edges_,
            faces_,
            cells_,
            owner_,
            scalField,
            lablField
        );
    }
}


// Return the status report interval
scalar dynamicTopoFvMesh::reportInterval() const
{
    // Default to 1 second
    scalar interval = 1.0;

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
    Pout<< "Checking index ranges...";

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
            Pout<< "Edge " << edgeI
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
            Pout<< "Edge " << edgeI
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
            // This might be a deleted face
            if (faceI < nOldFaces_)
            {
                if (reverseFaceMap_[faceI] == -1)
                {
                    continue;
                }
            }
            else
            if (deletedFaces_.found(faceI))
            {
                continue;
            }

            Pout<< "Face " << faceI
                << " has no vertex labels."
                << endl;

            nFailedChecks++;
            continue;
        }

        if (min(curFace) < 0 || max(curFace) > (points_.size()-1))
        {
            Pout<< "Face " << faceI
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
                Pout<< "Face " << faceI
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
                Pout<< "Cells for face: " << faceI << "::" << curFace << nl
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
            Pout<< "Cell " << cellI
                << " contains face labels out of range: " << curCell
                << " Max face index = " << (faces_.size()-1) << endl;

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
                Pout<< "Cell " << cellI
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

    Pout<< "Done." << endl;

    Pout<< "Checking face-cell connectivity...";

    forAll(nCellsPerFace, faceI)
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
        if (deletedFaces_.found(faceI))
        {
            continue;
        }

        // Determine patch
        label uPatch = whichPatch(faceI);

        if (nCellsPerFace[faceI] == 0)
        {
            // Looks like this is really an unused face.
            Pout<< "Face " << faceI << " :: " << faces_[faceI]
                << " is unused. Patch: "
                << (uPatch > -1 ? boundaryMesh()[uPatch].name() : "Internal")
                << nl << endl;

            nFailedChecks++;

            ConnectivityWarning()
                << nl << "Cell-Face connectivity is inconsistent."
                << endl;
        }
        else
        if (nCellsPerFace[faceI] != 2 && uPatch == -1)
        {
            // Internal face is not shared by exactly two cells
            Pout<< "Internal Face " << faceI
                << " :: " << faces_[faceI]
                << " Owner: " << owner_[faceI]
                << " Neighbour: " << neighbour_[faceI]
                << " is multiply connected." << nl
                << " nCellsPerFace: " << nCellsPerFace[faceI] << nl
                << " Patch: Internal" << nl
                << endl;

            // Loop through cells and find another instance
            forAll(cells_, cellI)
            {
                if (findIndex(cells_[cellI], faceI) > -1)
                {
                    Pout<< "  Cell: " << cellI
                        << "  :: " << cells_[cellI]
                        << endl;
                }
            }

            nFailedChecks++;

            ConnectivityWarning()
                << nl << "Cell-Face connectivity is inconsistent."
                << endl;
        }
        else
        if (nCellsPerFace[faceI] != 1 && uPatch > -1)
        {
            // Boundary face is not shared by exactly one cell
            Pout<< "Boundary Face " << faceI
                << " :: " << faces_[faceI]
                << " Owner: " << owner_[faceI]
                << " Neighbour: " << neighbour_[faceI]
                << " is multiply connected." << nl
                << " nCellsPerFace: " << nCellsPerFace[faceI] << nl
                << " Patch: " << boundaryMesh()[uPatch].name() << nl
                << endl;

            // Loop through cells and find another instance
            forAll(cells_, cellI)
            {
                if (findIndex(cells_[cellI], faceI) > -1)
                {
                    Pout<< "  Cell: " << cellI
                        << "  :: " << cells_[cellI]
                        << endl;
                }
            }

            nFailedChecks++;

            ConnectivityWarning()
                << nl << "Cell-Face connectivity is inconsistent."
                << endl;
        }
    }

    Pout<< "Done." << endl;

    Pout<< "Checking for unused points...";

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
            Pout<< nl << nl << "Point " << pointI
                << " is unused. " << endl;

            nFailedChecks++;

            ConnectivityWarning()
                << nl << "Point-Face connectivity is inconsistent."
                << endl;
        }
    }

    Pout<< "Done." << endl;

    Pout<< "Checking edge-face connectivity...";

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
                Pout<< nl << nl << "Edge: " << faceEdges[edgeI]
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
    dynamicLabelList bPatchIDs(10);
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
            Pout<< nl << nl << "Edge: " << edgeI << ": " << edges_[edgeI]
                << ": edgeFaces: " << edgeFaces
                << nl << UIndirectList<face>(faces_, edgeFaces)
                << nl << " Expected nFaces: " << nEdgeFaces[edgeI]
                << endl;

            nFailedChecks++;

            ConnectivityWarning()
                << nl << "Edge-Face connectivity is inconsistent."
                << endl;
        }

        label nBF = 0;
        bPatchIDs.clear();

        // Check if this edge belongs to faceEdges for each face
        forAll(edgeFaces, faceI)
        {
            const labelList& faceEdges = faceEdges_[edgeFaces[faceI]];

            if (findIndex(faceEdges, edgeI) == -1)
            {
                Pout<< nl << nl << "Edge: " << edgeI << ": " << edges_[edgeI]
                    << ", edgeFaces: " << edgeFaces
                    << nl << UIndirectList<face>(faces_, edgeFaces)
                    << "was not found in faceEdges of face: "
                    << edgeFaces[faceI] << ": " << faces_[edgeFaces[faceI]]
                    << nl << "faceEdges: " << faceEdges
                    << nl << UIndirectList<edge>(edges_, faceEdges)
                    << endl;

                nFailedChecks++;

                ConnectivityWarning()
                    << nl << "Edge-Face connectivity is inconsistent."
                    << endl;
            }

            if (neighbour_[edgeFaces[faceI]] == -1)
            {
                // Add to list of patch IDs
                bPatchIDs.append(whichPatch(edgeFaces[faceI]));

                nBF++;
            }
        }

        if (nBF == 0)
        {
            nInternalEdges++;

            // Check if this edge is actually internal.
            if (whichEdgePatch(edgeI) >= 0)
            {
                Pout<< "Edge: " << edgeI
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
                Pout<< "Edge: " << edgeI
                    << ": " << edges_[edgeI]
                    << " is on a boundary, but patch is specified as: "
                    << patchID << endl;

                nFailedChecks++;
            }
            else
            {
                patchInfo[patchID]++;
            }

            bool failedManifoldCheck = false;

            if (nBF > 2)
            {
                if (Pstream::parRun())
                {
                    // Pinched manifolds should be allowed in parallel
                    failedManifoldCheck = false;
                }
                else
                {
                    failedManifoldCheck = true;
                }
            }

            if (failedManifoldCheck)
            {
                // Write out for post-processing
                forAll(bPatchIDs, faceI)
                {
                    if (bPatchIDs[faceI] == -1)
                    {
                        Pout<< " Edge: " << edgeI
                            << " Face Patch: Internal" << nl;
                    }
                    else
                    {
                        Pout<< " Edge: " << edgeI
                            << " Face Patch: "
                            << boundaryMesh()[bPatchIDs[faceI]].name() << nl;
                    }
                }

                Pout<< endl;

                writeVTK("pinched_" + Foam::name(edgeI), edgeFaces, 2);

                Pout<< "Edge: " << edgeI
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
        Pout<< nl << "Internal edge-count is inconsistent." << nl
            << " Counted internal edges: " << nInternalEdges
            << " Actual count: " << nInternalEdges_ << endl;

        nFailedChecks++;
    }

    forAll(patchInfo, patchI)
    {
        if (patchInfo[patchI] != edgePatchSizes_[patchI])
        {
            Pout<< "Patch-count is inconsistent." << nl
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
            Pout<< nl << nl << "Edge: " << key
                << "::" << edges_[key]
                << ", edgeFaces: " << edgeFaces
                << " is internal, but contains boundary faces."
                << endl;

            nFailedChecks++;
        }

        if ((patch >= 0) && (nBF != 2))
        {
            if (!Pstream::parRun())
            {
                Pout<< nl << nl << "Edge: " << key
                    << "::" << edges_[key]
                    << ", edgeFaces: " << edgeFaces
                    << " is on a boundary patch, but doesn't contain"
                    << " two boundary faces."
                    << endl;

                nFailedChecks++;
            }
        }
    }

    Pout<< "Done." << endl;

    // Check coupled-patch sizes
    forAll(patchCoupling_, patchI)
    {
        if (patchCoupling_(patchI))
        {
            const coupleMap& cMap = patchCoupling_[patchI].map();

            label mSize = patchSizes_[cMap.masterIndex()];
            label sSize = patchSizes_[cMap.slaveIndex()];

            if (mSize != sSize)
            {
                Pout<< "Coupled patch-count is inconsistent." << nl
                    << " Master Patch: " << cMap.masterIndex()
                    << " Count: " << mSize << nl
                    << " Slave Patch: " << cMap.slaveIndex()
                    << " Count: " << sSize
                    << endl;

                nFailedChecks++;
            }
        }
    }

    if (is3D())
    {
        Pout<< "Checking point-edge connectivity...";

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
                    Pout<< nl << nl << "Point: " << pointI << nl
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
                Pout<< nl << nl << "Point: " << pointI << nl
                    << "pointEdges: " << pointEdges << nl
                    << "hlPointEdges: " << hlPointEdges[pointI]
                    << endl;

                nFailedChecks++;

                ConnectivityWarning()
                    << nl << "Size inconsistency."
                    << nl << "Point-Edge connectivity is inconsistent."
                    << endl;
            }
        }

        Pout<< "Done." << endl;
    }

    Pout<< "Checking cell-point connectivity...";

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
        // Check for hexahedral cells
        if
        (
            (cellToNode[cellI].size() == 8) &&
            (cells_[cellIndex[cellI]].size() == 6)
        )
        {
            continue;
        }

        if
        (
            (cellToNode[cellI].size() != 6 && is2D()) ||
            (cellToNode[cellI].size() != 4 && is3D())
        )
        {
            Pout<< nl << "Warning: Cell: "
                << cellIndex[cellI] << " is inconsistent. "
                << endl;

            const cell& failedCell = cells_[cellIndex[cellI]];

            Pout<< "Cell faces: " << failedCell << endl;

            forAll(failedCell, faceI)
            {
                Pout<< "\tFace: " << failedCell[faceI]
                    << " :: " << faces_[failedCell[faceI]]
                    << endl;

                const labelList& fEdges = faceEdges_[failedCell[faceI]];

                forAll(fEdges, edgeI)
                {
                    Pout<< "\t\tEdge: " << fEdges[edgeI]
                        << " :: " << edges_[fEdges[edgeI]]
                        << endl;
                }
            }

            nFailedChecks++;
        }
    }

    Pout<< "Done." << endl;

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


// Utility method to check whether the faces in triFaces will yield
// valid triangles when 'pointIndex' is moved to 'newPoint'.
//  - The routine performs metric-based checks.
//  - Returns 'true' if the collapse in NOT feasible.
//  - Does not reference member data, because this function
//    is also used on subMeshes
bool dynamicTopoFvMesh::checkCollapse
(
    const dynamicTopoFvMesh& mesh,
    const labelList& triFaces,
    const FixedList<label,2>& c0BdyIndex,
    const FixedList<label,2>& c1BdyIndex,
    const FixedList<label,2>& pointIndex,
    const FixedList<point,2>& newPoint,
    const FixedList<point,2>& oldPoint,
    scalar& collapseQuality,
    const bool checkNeighbour,
    bool forceOp
)
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

        const face& checkFace = mesh.faces_[triFaces[indexI]];

        // Configure a triangle face
        FixedList<point, 3> tFNew(vector::zero);
        FixedList<point, 3> tFOld(vector::zero);

        // Make necessary replacements
        forAll(checkFace, pointI)
        {
            tFNew[pointI] = mesh.points_[checkFace[pointI]];
            tFOld[pointI] = mesh.oldPoints_[checkFace[pointI]];

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
    if (collapseQuality < mesh.sliverThreshold_ && !forceOp)
    {
        if (debug > 3)
        {
            Pout<< " * * * 2D checkCollapse * * * " << nl
                << " collapseQuality: " << collapseQuality
                << " below threshold: " << mesh.sliverThreshold_
                << endl;
        }

        return true;
    }

    // Negative quality is a no-no
    if (collapseQuality < 0.0)
    {
        if (forceOp)
        {
            Pout<< " * * * 2D checkCollapse * * * " << nl
                << " Negative collapseQuality: " << collapseQuality << nl
                << " Operation cannot be forced."
                << endl;
        }

        return true;
    }

    // Negative old-area is also a no-no
    if (minArea < 0.0)
    {
        if (forceOp)
        {
            Pout<< " * * * 2D checkCollapse * * * " << nl
                << " minArea: " << minArea << nl
                << " Operation cannot be forced."
                << endl;
        }

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
    dynamicLabelList& cellsChecked,
    scalar& collapseQuality,
    bool forceOp
) const
{
    label faceIndex = -1;
    scalar cQuality = 0.0, oldVolume = 0.0, newVolume = 0.0;
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

        newVolume =
        (
            tetPointRef
            (
                points_[faceToCheck[2]],
                points_[faceToCheck[1]],
                points_[faceToCheck[0]],
                newPoint
            ).mag()
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

        newVolume =
        (
            tetPointRef
            (
                points_[faceToCheck[0]],
                points_[faceToCheck[1]],
                points_[faceToCheck[2]],
                newPoint
            ).mag()
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
        if (debug > 4)
        {
            Pout<< " * * * 3D checkCollapse * * * " << nl
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
            Pout<< " * * * 3D checkCollapse * * * " << nl
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
    if (oldVolume < 0.0 || (mag(oldVolume) < mag(0.1 * newVolume)))
    {
        if (forceOp)
        {
            Pout<< " * * * 3D checkCollapse * * * " << nl
                << "\nCollapsing cell: " << cellIndex
                << " containing points:\n"
                << faceToCheck[0] << "," << faceToCheck[1] << ","
                << faceToCheck[2] << "," << pointIndex << nl
                << "will yield an old-volume of: " << oldVolume
                << ", when " << pointIndex
                << " is moved to location: " << nl
                << oldPoint << nl
                << "newVolume: " << newVolume << nl
                << "Operation cannot be forced."
                << endl;
        }

        return true;
    }

    // No problems, so a collapse is feasible
    cellsChecked.append(cellIndex);

    // Update input quality
    collapseQuality = Foam::min(collapseQuality, cQuality);

    return false;
}


} // End namespace Foam

// ************************************************************************* //
