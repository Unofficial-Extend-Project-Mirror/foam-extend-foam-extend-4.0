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
    An implementation of dynamic changes to mesh-topology

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "dynamicTopoFvMesh.H"
#include "addToRunTimeSelectionTable.H"

#include "eMesh.H"
#include "Stack.H"
#include "triFace.H"
#include "changeMap.H"
#include "clockTime.H"
#include "volFields.H"
#include "MeshObject.H"
#include "topoMapper.H"
#include "mapPolyMesh.H"
#include "MapFvFields.H"
#include "SortableList.H"
#include "motionSolver.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"
#include "lengthScaleEstimator.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynamicTopoFvMesh,0);

addToRunTimeSelectionTable(dynamicFvMesh, dynamicTopoFvMesh, IOobject);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from IOobject
dynamicTopoFvMesh::dynamicTopoFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    topoChangeFlag_(false),
    dict_
    (
        IOobject
        (
            "dynamicMeshDict",
            polyMesh::time().constant(),
            (*this),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mandatory_
    (
        dict_.subDict("dynamicTopoFvMesh").lookup("allOptionsMandatory")
    ),
    twoDMesh_(polyMesh::nGeometricD() == 2 ? true : false),
    edgeRefinement_
    (
        dict_.subDict("dynamicTopoFvMesh").lookup("edgeRefinement")
    ),
    bandWidthReduction_(false),
    interval_(1),
    eMeshPtr_(NULL),
    mapper_(NULL),
    motionSolver_(NULL),
    lengthEstimator_(NULL),
    oldPoints_(polyMesh::points()),
    points_(polyMesh::points()),
    faces_(polyMesh::faces()),
    owner_(polyMesh::faceOwner()),
    neighbour_(polyMesh::faceNeighbour()),
    cells_(primitiveMesh::cells()),
    oldPatchSizes_(boundaryMesh().size(),0),
    patchSizes_(boundaryMesh().size(),0),
    oldPatchStarts_(boundaryMesh().size(),-1),
    patchStarts_(boundaryMesh().size(),-1),
    oldEdgePatchSizes_(boundaryMesh().size(),0),
    edgePatchSizes_(boundaryMesh().size(),0),
    oldEdgePatchStarts_(boundaryMesh().size(),-1),
    edgePatchStarts_(boundaryMesh().size(),-1),
    oldPatchNMeshPoints_(boundaryMesh().size(),-1),
    patchNMeshPoints_(boundaryMesh().size(),-1),
    nOldPoints_(primitiveMesh::nPoints()),
    nPoints_(primitiveMesh::nPoints()),
    nOldEdges_(0),
    nEdges_(0),
    nOldFaces_(primitiveMesh::nFaces()),
    nFaces_(primitiveMesh::nFaces()),
    nOldCells_(primitiveMesh::nCells()),
    nCells_(primitiveMesh::nCells()),
    nOldInternalFaces_(primitiveMesh::nInternalFaces()),
    nInternalFaces_(primitiveMesh::nInternalFaces()),
    nOldInternalEdges_(0),
    nInternalEdges_(0),
    maxModifications_(-1),
    statistics_(0),
    sliverThreshold_(0.1),
    slicePairs_(0),
    maxTetsPerEdge_(-1),
    swapDeviation_(0.0),
    allowTableResize_(false)
{
    // Check the size of owner/neighbour
    if (owner_.size() != neighbour_.size())
    {
        // Size up to number of faces
        neighbour_.setSize(nFaces_, -1);
    }

    const polyBoundaryMesh& boundary = polyMesh::boundaryMesh();

    // Initialize the multiThreading environment
    initializeThreadingEnvironment();

    // Read optional parameters.
    readOptionalParameters();

    // Initialize patch-size information
    for (label i = 0; i < boundaryMesh().size(); i++)
    {
        patchNMeshPoints_[i] = boundary[i].meshPoints().size();
        oldPatchSizes_[i] = patchSizes_[i]  = boundary[i].size();
        oldPatchStarts_[i] = patchStarts_[i] = boundary[i].start();
        oldPatchNMeshPoints_[i] = patchNMeshPoints_[i];
    }

    // Open the tetMetric dynamic-link library (for 3D only)
    loadMetric();

    // Initialize edge-related connectivity structures
    initEdges();

    // Load the mesh-motion solver
    loadMotionSolver();

    // Load the field-mapper
    loadFieldMapper();

    // Set sizes for the reverse maps
    reversePointMap_.setSize(nPoints_, -7);
    reverseEdgeMap_.setSize(nEdges_, -7);
    reverseFaceMap_.setSize(nFaces_, -7);
    reverseCellMap_.setSize(nCells_, -7);

    // Load the length-scale estimator,
    // and read refinement options
    loadLengthScaleEstimator();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dynamicTopoFvMesh::~dynamicTopoFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Insert the specified cell to the mesh.
label dynamicTopoFvMesh::insertCell
(
    const cell& newCell,
    const scalar lengthScale,
    const label zoneID
)
{
    label newCellIndex = cells_.size();

    if (debug > 2)
    {
        Info << "Inserting cell: "
             << newCellIndex << ": "
             << newCell << endl;
    }

    cells_.append(newCell);

    if (edgeRefinement_)
    {
        lengthScale_.append(lengthScale);
    }

    // Add to the zone if necessary
    if (zoneID >= 0)
    {
        addedCellZones_.insert(newCellIndex, zoneID);
    }

    nCells_++;

    return newCellIndex;
}


// Remove the specified cell from the mesh
void dynamicTopoFvMesh::removeCell
(
    const label cIndex
)
{
    if (debug > 2)
    {
        Info << "Removing cell: "
             << cIndex << ": "
             << cells_[cIndex]
             << endl;
    }

    cells_[cIndex].clear();

    if (edgeRefinement_)
    {
        lengthScale_[cIndex] = -1.0;
    }

    // Update the number of cells, and the reverse cell map
    nCells_--;

    if (cIndex < nOldCells_)
    {
        reverseCellMap_[cIndex] = -1;
    }
    else
    {
        // Store this information for the reOrdering stage
        deletedCells_.insert(cIndex);
    }

    // Check if this cell was added to a zone
    if (addedCellZones_.found(cIndex))
    {
        addedCellZones_.erase(cIndex);
    }

    // Check if the cell was added in the current morph, and delete
    forAll(cellsFromPoints_, indexI)
    {
        if (cellsFromPoints_[indexI].index() == cIndex)
        {
            // Remove entry from the list
            meshOps::removeIndex(indexI, cellsFromPoints_);
            break;
        }
    }

    forAll(cellsFromEdges_, indexI)
    {
        if (cellsFromEdges_[indexI].index() == cIndex)
        {
            // Remove entry from the list
            meshOps::removeIndex(indexI, cellsFromEdges_);
            break;
        }
    }

    forAll(cellsFromFaces_, indexI)
    {
        if (cellsFromFaces_[indexI].index() == cIndex)
        {
            // Remove entry from the list
            meshOps::removeIndex(indexI, cellsFromFaces_);
            break;
        }
    }

    forAll(cellsFromCells_, indexI)
    {
        if (cellsFromCells_[indexI].index() == cIndex)
        {
            // Remove entry from the list
            meshOps::removeIndex(indexI, cellsFromCells_);
            break;
        }
    }
}


// Utility method for face-insertion
label dynamicTopoFvMesh::insertFace
(
    const label patch,
    const face& newFace,
    const label newOwner,
    const label newNeighbour,
    const label zoneID
)
{
    // Append the specified face to each face-related list.
    // Reordering is performed after all pending changes
    // (flips, bisections, contractions, etc) have been made to the mesh
    label newFaceIndex = faces_.size();

    faces_.append(newFace);
    owner_.append(newOwner);
    neighbour_.append(newNeighbour);

    if (debug > 2)
    {
        Info << "Inserting face: "
             << newFaceIndex << ": "
             << newFace
             << " Owner: " << newOwner
             << " Neighbour: " << newNeighbour;

        Info << " Patch: ";

        if (patch == -1)
        {
            Info << "Internal" << endl;
        }
        else
        {
            Info << boundaryMesh()[patch].name() << endl;
        }
    }

    // Keep track of added boundary faces in a separate hash-table
    // This information will be required at the reordering stage
    addedFacePatches_.insert(newFaceIndex,patch);

    if (newNeighbour == -1)
    {
        // Modify patch information for this boundary face
        patchSizes_[patch]++;

        for (label i = (patch + 1); i < boundaryMesh().size(); i++)
        {
            patchStarts_[i]++;
        }
    }
    else
    {
        // Increment the number of internal faces,
        // and subsequent patch-starts
        nInternalFaces_++;

        for (label i = 0; i < boundaryMesh().size(); i++)
        {
            patchStarts_[i]++;
        }
    }

    // Add to the zone if explicitly specified
    if (zoneID >= 0)
    {
        addedFaceZones_.insert(newFaceIndex, zoneID);
    }

    // Increment the total face count
    nFaces_++;

    return newFaceIndex;
}


// Remove the specified face from the mesh
void dynamicTopoFvMesh::removeFace
(
    const label fIndex
)
{
    if (debug > 2)
    {
        Info << "Removed face: "
             << fIndex << ": "
             << faces_[fIndex] << endl;
    }

    // Identify the patch for this face
    label patch = whichPatch(fIndex);

    if (patch >= 0)
    {
        // Modify patch information for this boundary face
        patchSizes_[patch]--;

        for (label i = (patch + 1); i < boundaryMesh().size(); i++)
        {
            patchStarts_[i]--;
        }
    }
    else
    {
        // Decrement the internal face count, and subsequent patch-starts
        nInternalFaces_--;

        forAll(patchStarts_, patchI)
        {
            patchStarts_[patchI]--;
        }
    }

    // Clear entities.
    faces_[fIndex].clear();
    owner_[fIndex] = -1;
    neighbour_[fIndex] = -1;
    faceEdges_[fIndex].clear();

    if (twoDMesh_)
    {
        // Remove from the stack as well
        forAll(entityStack_, stackI)
        {
            stack(stackI).remove(fIndex);
        }
    }

    // Update the reverse face-map, but only if this is a face that existed
    // at time [n]. Added faces which are deleted during the topology change
    // needn't be updated.
    if (fIndex < nOldFaces_)
    {
        reverseFaceMap_[fIndex] = -1;
    }
    else
    {
        // Store this information for the reOrdering stage
        deletedFaces_.insert(fIndex);
    }

    // Check and remove from the list of added face patches
    if (addedFacePatches_.found(fIndex))
    {
        addedFacePatches_.erase(fIndex);
    }

    // Check if this face was added to a zone
    if (addedFaceZones_.found(fIndex))
    {
        addedFaceZones_.erase(fIndex);
    }

    // Check if the face was added in the current morph, and delete
    forAll(facesFromPoints_, indexI)
    {
        if (facesFromPoints_[indexI].index() == fIndex)
        {
            // Remove entry from the list
            meshOps::removeIndex
            (
                indexI,
                facesFromPoints_
            );

            break;
        }
    }

    forAll(facesFromEdges_, indexI)
    {
        if (facesFromEdges_[indexI].index() == fIndex)
        {
            // Remove entry from the list
            meshOps::removeIndex
            (
                indexI,
                facesFromEdges_
            );

            break;
        }
    }

    forAll(facesFromFaces_, indexI)
    {
        if (facesFromFaces_[indexI].index() == fIndex)
        {
            // Remove entry from the list
            meshOps::removeIndex
            (
                indexI,
                facesFromFaces_
            );

            break;
        }
    }

    // Remove from the flipFaces list, if necessary
    if (flipFaces_.found(fIndex))
    {
        flipFaces_.erase(fIndex);
    }

    // Decrement the total face-count
    nFaces_--;
}


// Insert the specified edge to the mesh
label dynamicTopoFvMesh::insertEdge
(
    const label patch,
    const edge& newEdge,
    const labelList& edgeFaces,
    const labelList& edgePoints
)
{
    label newEdgeIndex = edges_.size();

    edges_.append(newEdge);
    edgeFaces_.append(edgeFaces);

    if (!twoDMesh_)
    {
        edgePoints_.append(edgePoints);
    }

    if (debug > 2)
    {
        Info << "Inserting edge: "
             << newEdgeIndex << ": "
             << newEdge;

        Info << " Patch: ";

        if (patch == -1)
        {
            Info << "Internal" << endl;
        }
        else
        {
            Info << boundaryMesh()[patch].name() << endl;
        }

        if (!twoDMesh_)
        {
            if (findIndex(edgePoints, -1) != -1)
            {
                FatalErrorIn
                (
                    "\n"
                    "label dynamicTopoFvMesh::insertEdge\n"
                    "(\n"
                    "    const label patch,\n"
                    "    const edge& newEdge,\n"
                    "    const labelList& edgeFaces,\n"
                    "    const labelList& edgePoints\n"
                    ")"
                )
                    << " EdgePoints is incorrectly specified." << nl
                    << " edgePoints: " << edgePoints << nl
                    << abort(FatalError);
            }
        }
    }

    // Keep track of added edges in a separate hash-table
    // This information will be required at the reordering stage
    addedEdgePatches_.insert(newEdgeIndex,patch);

    if (patch >= 0)
    {
        // Modify patch information for this boundary edge
        edgePatchSizes_[patch]++;

        for (label i = (patch + 1); i < boundaryMesh().size(); i++)
        {
            edgePatchStarts_[i]++;
        }
    }
    else
    {
        // Increment the number of internal edges, and subsequent patch-starts
        nInternalEdges_++;

        for (label i = 0; i < boundaryMesh().size(); i++)
        {
            edgePatchStarts_[i]++;
        }
    }

    // Size-up the pointEdges list
    if (!twoDMesh_)
    {
        meshOps::sizeUpList(newEdgeIndex, pointEdges_[newEdge[0]]);
        meshOps::sizeUpList(newEdgeIndex, pointEdges_[newEdge[1]]);
    }

    // Increment the total edge count
    nEdges_++;

    return newEdgeIndex;
}


// Remove the specified edge from the mesh
void dynamicTopoFvMesh::removeEdge
(
    const label eIndex
)
{
    if (debug > 2)
    {
        Info << "Removing edge: "
             << eIndex << ": "
             << edges_[eIndex] << endl;
    }

    if (!twoDMesh_)
    {
        // Remove the edgePoints entry
        edgePoints_[eIndex].clear();

        // Size-down the pointEdges list
        if (pointEdges_[edges_[eIndex][0]].size())
        {
            meshOps::sizeDownList(eIndex, pointEdges_[edges_[eIndex][0]]);
        }

        if (pointEdges_[edges_[eIndex][1]].size())
        {
            meshOps::sizeDownList(eIndex, pointEdges_[edges_[eIndex][1]]);
        }

        // Remove from the stack as well
        forAll(entityStack_, stackI)
        {
            stack(stackI).remove(eIndex);
        }
    }

    edges_[eIndex] = edge(-1, -1);
    edgeFaces_[eIndex].clear();

    // Identify the patch for this edge
    label patch = whichEdgePatch(eIndex);

    if (patch >= 0)
    {
        // Modify patch information for this boundary edge
        edgePatchSizes_[patch]--;

        for (label i = (patch + 1); i < boundaryMesh().size(); i++)
        {
            edgePatchStarts_[i]--;
        }
    }
    else
    {
        // Decrement the internal edge count, and subsequent patch-starts
        nInternalEdges_--;

        forAll(edgePatchStarts_, patchI)
        {
            edgePatchStarts_[patchI]--;
        }
    }

    // Update reverse edge-map, but only if this is an edge that existed
    // at time [n]. Added edges which are deleted during the topology change
    // needn't be updated.
    if (eIndex < nOldEdges_)
    {
        reverseEdgeMap_[eIndex] = -1;
    }
    else
    {
        // Store this information for the reOrdering stage
        deletedEdges_.insert(eIndex);
    }

    // Check and remove from the list of added edge patches
    if (addedEdgePatches_.found(eIndex))
    {
        addedEdgePatches_.erase(eIndex);
    }

    // Decrement the total edge-count
    nEdges_--;
}


// Insert the specified point to the mesh
label dynamicTopoFvMesh::insertPoint
(
    const point& newPoint,
    const point& oldPoint,
    const labelList& mapPoints,
    const label zoneID
)
{
    // Add a new point to the end of the list
    label newPointIndex = points_.size();

    points_.append(newPoint);
    oldPoints_.append(oldPoint);

    if (debug > 2)
    {
        Info << "Inserting new point: "
             << newPointIndex << ": "
             << newPoint
             << " and old point: "
             << oldPoint
             << "  Mapped from: "
             << mapPoints << endl;
    }

    // Make a pointsFromPoints entry
    meshOps::sizeUpList
    (
        objectMap(newPointIndex, mapPoints),
        pointsFromPoints_
    );

    // Add an empty entry to pointEdges as well.
    // This entry can be sized-up appropriately at a later stage.
    if (!twoDMesh_)
    {
        pointEdges_.append(labelList(0));
    }

    // Add to the zone if necessary
    if (zoneID >= 0)
    {
        addedPointZones_.insert(newPointIndex, zoneID);
    }

    nPoints_++;

    return newPointIndex;
}


// Remove the specified point from the mesh
void dynamicTopoFvMesh::removePoint
(
    const label pIndex
)
{
    if (debug > 2)
    {
        Info << "Removing point: "
             << pIndex << ": "
             << points_[pIndex] << endl;
    }

    // Remove the point
    // (or just make sure that it's never used anywhere else)
    // points_[pIndex] = point();
    // oldPoints_[pIndex] = point();

    // Remove pointEdges as well
    if (!twoDMesh_)
    {
        pointEdges_[pIndex].clear();
    }

    // Update the reverse point map
    if (pIndex < nOldPoints_)
    {
        reversePointMap_[pIndex] = -1;
    }
    else
    {
        deletedPoints_.insert(pIndex);
    }

    // Check if this point was added to a zone
    if (addedPointZones_.found(pIndex))
    {
        addedPointZones_.erase(pIndex);
    }

    // Check if the point was added in the current morph, and delete
    forAll(pointsFromPoints_, indexI)
    {
        if (pointsFromPoints_[indexI].index() == pIndex)
        {
            // Remove entry from the list
            meshOps::removeIndex(indexI, pointsFromPoints_);
            break;
        }
    }

    // Decrement the total point-count
    nPoints_--;
}


// Utility method to build edgePoints for an edge [3D].
// Assumes that edgeFaces information is consistent.
void dynamicTopoFvMesh::buildEdgePoints
(
    const label eIndex,
    const label checkIndex
)
{
    bool found = false;
    label faceIndex = -1, cellIndex = -1;
    label otherPoint = -1, nextPoint = -1;

    // Obtain references
    const edge& edgeToCheck = edges_[eIndex];
    const labelList& eFaces = edgeFaces_[eIndex];

    // Re-size the list first
    labelList& ePoints = edgePoints_[eIndex];
    ePoints.setSize(eFaces.size(), -1);

    if (whichEdgePatch(eIndex) == -1)
    {
        // Internal edge.
        // Pick the first face and start with that
        faceIndex = eFaces[0];
    }
    else
    {
        // Need to find a properly oriented start-face
        forAll(eFaces, faceI)
        {
            if (whichPatch(eFaces[faceI]) > -1)
            {
                meshOps::findIsolatedPoint
                (
                    faces_[eFaces[faceI]],
                    edgeToCheck,
                    otherPoint,
                    nextPoint
                );

                if (nextPoint == edgeToCheck[checkIndex])
                {
                    faceIndex = eFaces[faceI];
                    break;
                }
            }
        }
    }

    // Shuffle vertices to appear in CCW order
    forAll(ePoints, indexI)
    {
        meshOps::findIsolatedPoint
        (
            faces_[faceIndex],
            edgeToCheck,
            otherPoint,
            nextPoint
        );

        // Add the isolated point
        ePoints[indexI] = otherPoint;

        // Figure out how this edge is oriented.
        if (nextPoint == edgeToCheck[checkIndex])
        {
            // Counter-clockwise. Pick the owner.
            cellIndex = owner_[faceIndex];
        }
        else
        if (whichPatch(faceIndex) == -1)
        {
            // Clockwise. Pick the neighbour.
            cellIndex = neighbour_[faceIndex];
        }
        else
        {
            // Looks like we've hit a boundary face. Break out.
            break;
        }

        const cell& cellToCheck = cells_[cellIndex];

        found = false;

        // Assuming tet-cells,
        // Loop through edgeFaces and get the next face
        forAll(eFaces, faceI)
        {
            if
            (
                eFaces[faceI] != faceIndex
             && eFaces[faceI] == cellToCheck[0]
            )
            {
                faceIndex = cellToCheck[0];
                found = true; break;
            }

            if
            (
                eFaces[faceI] != faceIndex
             && eFaces[faceI] == cellToCheck[1]
            )
            {
                faceIndex = cellToCheck[1];
                found = true; break;
            }

            if
            (
                eFaces[faceI] != faceIndex
             && eFaces[faceI] == cellToCheck[2]
            )
            {
                faceIndex = cellToCheck[2];
                found = true; break;
            }

            if
            (
                eFaces[faceI] != faceIndex
             && eFaces[faceI] == cellToCheck[3]
            )
            {
                faceIndex = cellToCheck[3];
                found = true; break;
            }
        }

        if (!found)
        {
            Info << "edgeFaces: " << endl;

            forAll(eFaces, faceI)
            {
                Info << " Face: " << eFaces[faceI]
                     << ":: " << faces_[eFaces[faceI]]
                     << " Owner: " << owner_[eFaces[faceI]]
                     << " Neighbour: " << neighbour_[eFaces[faceI]]
                     << endl;
            }

            writeVTK("vRingEdgeFaces", eFaces, 2);
            writeVTK("vRingCellToCheck", cellIndex);

            // Something's terribly wrong
            FatalErrorIn
            (
                "\n"
                "void dynamicTopoFvMesh::buildEdgePoints\n"
                "(\n"
                "    const label eIndex,\n"
                "    const label checkIndex\n"
                ")"
            )
                << " Failed to determine a vertex ring. " << nl
                << " edgeFaces connectivity is inconsistent. " << nl
                << " Edge: " << eIndex << ":: " << edgeToCheck << nl
                << " edgeFaces: " << eFaces << nl
                << " Patch: " << whichEdgePatch(eIndex) << nl
                << " cellIndex: " << cellIndex
                << " :: " << cellToCheck << nl
                << " Current edgePoints: " << ePoints
                << abort(FatalError);
        }
    }
}


// Test an edge / face for proximity with other faces on proximity patches
// and return the scalar distance to an oppositely-oriented face.
scalar dynamicTopoFvMesh::testProximity
(
    const label index,
    label& proximityFace
) const
{
    scalar proxDistance = GREAT, testStep = 0.0;
    vector gCentre = vector::zero, gNormal = vector::zero;

    if (twoDMesh_)
    {
        // Obtain the face-normal.
        gNormal = faces_[index].normal(points_);

        // Obtain the face centre.
        gCentre = faces_[index].centre(points_);

        // Fetch the edge
        const edge& edgeToCheck = edges_[getTriBoundaryEdge(index)];

        // Calculate a test step-size
        testStep =
        (
            linePointRef
            (
                points_[edgeToCheck.start()],
                points_[edgeToCheck.end()]
            ).mag()
        );
    }
    else
    {
        const edge& edgeToCheck = edges_[index];
        const labelList& eFaces = edgeFaces_[index];

        linePointRef lpr
        (
            points_[edgeToCheck.start()],
            points_[edgeToCheck.end()]
        );

        // Obtain the edge centre.
        gCentre = lpr.centre();

        // Obtain the edge-normal
        forAll(eFaces, faceI)
        {
            if (neighbour_[eFaces[faceI]] == -1)
            {
                // Obtain the normal.
                gNormal += faces_[eFaces[faceI]].normal(points_);
            }
        }

        // Calculate a test step-size
        testStep = lpr.mag();
    }

    // Normalize
    gNormal /= (mag(gNormal) + VSMALL);

    // Test for proximity, and mark slice-pairs
    // is the distance is below threshold.
    if
    (
        lengthEstimator().testProximity
        (
            gCentre,
            gNormal,
            testStep,
            proximityFace,
            proxDistance
        )
    )
    {
        labelPair proxPoints(-1, -1);
        bool foundPoint = false;

        if (twoDMesh_)
        {
            proxPoints.first() = index;
            proxPoints.second() = proximityFace;

            if
            (
                (faces_[index].size() == 4) &&
                (polyMesh::faces()[proximityFace].size() == 4)
            )
            {
                foundPoint = true;
            }
        }
        else
        {
            const edge& thisEdge = edges_[index];
            const face& proxFace = polyMesh::faces()[proximityFace];

            // Check if any points on this face are still around.
            // If yes, mark one of them as the end point
            // for Dijkstra's algorithm. The start point will be a point
            // on this edge.
            proxPoints.first() = thisEdge[0];

            forAll(proxFace, pointI)
            {
                if (reversePointMap_[proxFace[pointI]] != -1)
                {
                    proxPoints.second() = proxFace[pointI];
                    foundPoint = true;
                    break;
                }
            }
        }

        if (foundPoint)
        {
            // Lock the edge mutex
            entityMutex_[1].lock();

            label newSize = slicePairs_.size() + 1;

            // Const-cast slicePairs for resize
            List<labelPair>& sP = const_cast<List<labelPair>&>(slicePairs_);

            // Add this entry as a candidate for mesh slicing.
            sP.setSize(newSize, proxPoints);

            // Unlock the edge mutex
            entityMutex_[1].unlock();
        }
    }

    return proxDistance;
}


// Calculate the edge length-scale for the mesh
void dynamicTopoFvMesh::calculateLengthScale(bool dump)
{
    if (!edgeRefinement_)
    {
        return;
    }

    Switch dumpLengthScale(false);

    const dictionary& meshDict = dict_.subDict("dynamicTopoFvMesh");

    if (meshDict.found("dumpLengthScale") || mandatory_)
    {
        dumpLengthScale = readBool(meshDict.lookup("dumpLengthScale"));
    }

    autoPtr<volScalarField> lsfPtr(NULL);

    if (dumpLengthScale && time().outputTime() && dump)
    {
        lsfPtr.set
        (
            new volScalarField
            (
                IOobject
                (
                    "lengthScale",
                    time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                *this,
                dimensionedScalar("scalar", dimLength, 0)
            )
        );
    }

    // Bail-out if a dumping was not requested in dictionary.
    if (dump && !dumpLengthScale)
    {
        return;
    }

    // Size the field and calculate length-scale
    lengthScale_.setSize(nCells_, 0.0);

    lengthEstimator().calculateLengthScale(lengthScale_);

    // Check if length-scale is to be dumped to disk.
    if (dumpLengthScale && time().outputTime() && dump)
    {
        // Obtain length-scale values from the mesh
        lsfPtr->internalField() = lengthScale_;

        lsfPtr->write();
    }
}


// Read optional dictionary parameters
void dynamicTopoFvMesh::readOptionalParameters(bool reRead)
{
    // Read from disk
    dict_.readIfModified();

    const dictionary& meshSubDict = dict_.subDict("dynamicTopoFvMesh");

    // Enable/disable run-time debug level
    if (meshSubDict.found("debug") || mandatory_)
    {
        debug = readLabel(meshSubDict.lookup("debug"));
    }
    else
    {
        debug = 0;
    }

    // Set debug option for underlying classes as well.
    lengthScaleEstimator::debug = debug;

    if (debug > 2)
    {
        fvMesh::debug = true;
        polyMesh::debug = true;
    }

    // Re-read edge-refinement options, if necessary
    if (edgeRefinement_ && reRead)
    {
        lengthEstimator().readRefinementOptions
        (
            meshSubDict.subDict("refinementOptions"),
            true,
            mandatory_
        );
    }

    // Read re-mesh interval
    if (meshSubDict.found("interval") || mandatory_)
    {
        interval_ = readLabel(meshSubDict.lookup("interval"));
    }
    else
    {
        interval_ = 1;
    }

    // Update bandwidth reduction switch
    if (meshSubDict.found("bandwidthReduction") || mandatory_)
    {
        bandWidthReduction_.readIfPresent("bandwidthReduction", meshSubDict);
    }

    // Update threshold for sliver cells
    if (meshSubDict.found("sliverThreshold") || mandatory_)
    {
        sliverThreshold_ = readScalar(meshSubDict.lookup("sliverThreshold"));

        if (sliverThreshold_ > 1.0 || sliverThreshold_ < 0.0)
        {
            FatalErrorIn("void dynamicTopoFvMesh::readOptionalParameters()")
                << " Sliver threshold out of range [0..1]"
                << abort(FatalError);
        }
    }

    // Update limit for max number of bisections / collapses
    if (meshSubDict.found("maxModifications") || mandatory_)
    {
        maxModifications_ = readLabel(meshSubDict.lookup("maxModifications"));
    }

    // Update limit for swap on curved surfaces
    if (meshSubDict.found("swapDeviation") || mandatory_)
    {
        swapDeviation_ = readScalar(meshSubDict.lookup("swapDeviation"));

        if (swapDeviation_ > 1.0 || swapDeviation_ < 0.0)
        {
            FatalErrorIn("void dynamicTopoFvMesh::readOptionalParameters()")
                << " Swap deviation out of range [0..1]"
                << abort(FatalError);
        }
    }

    // For tetrahedral meshes...
    if (!twoDMesh_)
    {
        // Check if swapping is to be avoided on any patches
        if (meshSubDict.found("noSwapPatches") || mandatory_)
        {
            wordList noSwapPatches =
            (
                meshSubDict.subDict("noSwapPatches").toc()
            );

            noSwapPatchIDs_.setSize(noSwapPatches.size());

            label indexI = 0;

            forAll(noSwapPatches, wordI)
            {
                const word& patchName = noSwapPatches[wordI];

                noSwapPatchIDs_[indexI++] =
                (
                    boundaryMesh().findPatchID(patchName)
                );
            }
        }

        // Check if a limit has been imposed on maxTetsPerEdge
        if (meshSubDict.found("maxTetsPerEdge") || mandatory_)
        {
            maxTetsPerEdge_ = readLabel(meshSubDict.lookup("maxTetsPerEdge"));
        }
        else
        {
            maxTetsPerEdge_ = 7;
        }

        // Check if programming tables can be resized at runtime
        if (meshSubDict.found("allowTableResize") || mandatory_)
        {
            allowTableResize_ =
            (
                readBool(meshSubDict.lookup("allowTableResize"))
            );
        }
        else
        {
            allowTableResize_ = false;
        }
    }
}


// Initialize edge related connectivity lists
void dynamicTopoFvMesh::initEdges()
{
    // Initialize eMesh, and copy to local lists
    eMeshPtr_.set(new eMesh(*this));

    // Obtain information
    nEdges_ = eMeshPtr_->nEdges();
    nInternalEdges_ = eMeshPtr_->nInternalEdges();
    edgePatchSizes_ = eMeshPtr_->boundary().patchSizes();
    edgePatchStarts_ = eMeshPtr_->boundary().patchStarts();

    // Set old edge information
    nOldEdges_ = nEdges_;
    nOldInternalEdges_ = nInternalEdges_;
    oldEdgePatchSizes_ = edgePatchSizes_;
    oldEdgePatchStarts_ = edgePatchStarts_;

    // Set local lists with edge connectivity information
    edges_ = eMeshPtr_->edges();
    edgeFaces_ = eMeshPtr_->edgeFaces();
    faceEdges_ = eMeshPtr_->faceEdges();

    if (!twoDMesh_)
    {
        pointEdges_ = eMeshPtr_->pointEdges();
        edgePoints_ = eMeshPtr_->edgePoints();
    }

    // Clear out unwanted eMesh connectivity
    eMeshPtr_->clearOut();
}


// Load the mesh-quality metric from the library
void dynamicTopoFvMesh::loadMetric()
{
    if (twoDMesh_)
    {
        return;
    }

    // Specify the dictionary we would be looking in...
    const dictionary& meshDict = dict_.subDict("dynamicTopoFvMesh");

    // Select an appropriate metric
    tetMetric_ = tetMetric::New(meshDict, meshDict.lookup("tetMetric"));
}


// Load the mesh-motion solver
void dynamicTopoFvMesh::loadMotionSolver()
{
    if (motionSolver_.valid())
    {
        FatalErrorIn
        (
            "void dynamicTopoFvMesh::loadMotionSolver() "
        ) << nl << " Motion solver already loaded. "
          << abort(FatalError);
    }
    else
    if (dict_.found("solver"))
    {
        motionSolver_ = motionSolver::New(*this);
    }
}


// Load the field mapper
void dynamicTopoFvMesh::loadFieldMapper()
{
    if (mapper_.valid())
    {
        FatalErrorIn
        (
            "void dynamicTopoFvMesh::loadFieldMapper() "
        ) << nl << " Field mapper already loaded. "
          << abort(FatalError);
    }
    else
    {
        mapper_.set(new topoMapper(*this, dict_));
    }
}


// Load the length scale estimator
void dynamicTopoFvMesh::loadLengthScaleEstimator()
{
    if (edgeRefinement_)
    {
        if (lengthEstimator_.valid())
        {
            FatalErrorIn
            (
                "void dynamicTopoFvMesh::loadLengthScaleEstimator() "
            ) << nl << " Length scale estimator already loaded. "
              << abort(FatalError);
        }
        else
        {
            lengthEstimator_.set
            (
                new lengthScaleEstimator(*this)
            );
        }

        // Read options
        lengthEstimator().readRefinementOptions
        (
            dict_.subDict("dynamicTopoFvMesh").subDict("refinementOptions"),
            false,
            mandatory_
        );
    }
}


// Initialize the threading environment.
//  - Provides an override option to avoid reading from the dictionary.
void dynamicTopoFvMesh::initializeThreadingEnvironment
(
    const label specThreads
)
{
    // Initialize an IOobject for the IOmultiThreader object
    IOobject io
    (
        "threader",
        this->time().timeName(),
        (*this),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    );

    if (specThreads > 0)
    {
        threader_.set(new IOmultiThreader(io, specThreads));
    }
    else
    {
        if (dict_.subDict("dynamicTopoFvMesh").found("threads") || mandatory_)
        {
            threader_.set
            (
                new IOmultiThreader
                (
                    io,
                    readLabel
                    (
                        dict_.subDict("dynamicTopoFvMesh").lookup("threads")
                    )
                )
            );
        }
        else
        {
            threader_.set(new IOmultiThreader(io, 1));
        }
    }

    // Get the number of threads and allocate threadHandlers
    label nThreads = threader_->getNumThreads();

    if (nThreads == 1)
    {
        handlerPtr_.setSize(1);

        handlerPtr_.set
        (
            0,
            new meshHandler(*this, threader())
        );

        handlerPtr_[0].setMaster();

        // Size the stacks
        entityStack_.setSize(1);
    }
    else
    {
        // Index '0' is master, rest are slaves
        handlerPtr_.setSize(nThreads + 1);

        // Size the stacks
        entityStack_.setSize(nThreads + 1);

        forAll(handlerPtr_, threadI)
        {
            handlerPtr_.set
            (
                threadI,
                new meshHandler(*this, threader())
            );

            if (threadI == 0)
            {
                handlerPtr_[0].setMaster();
            }
            else
            {
                handlerPtr_[threadI].setID(threader_->getID(threadI - 1));
                handlerPtr_[threadI].setSlave();
            }
        }
    }
}


// 2D Edge-swapping engine
void dynamicTopoFvMesh::swap2DEdges(void *argument)
{
    // Recast the argument
    meshHandler *thread = static_cast<meshHandler*>(argument);

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::START);
    }

    dynamicTopoFvMesh& mesh = thread->reference();

    // Figure out which thread this is...
    label tIndex = mesh.self();

    // Set the timer
    clockTime sTimer;

    bool reported = false;
    label stackSize = mesh.stack(tIndex).size();
    scalar interval = mesh.reportInterval(), oIndex = 0.0, nIndex = 0.0;

    oIndex = ::floor(sTimer.elapsedTime() / interval);

    // Pick items off the stack
    while (!mesh.stack(tIndex).empty())
    {
        // Report progress
        if (thread->master())
        {
            // Update the index, if its changed
            nIndex = ::floor(sTimer.elapsedTime() / interval);

            if ((nIndex - oIndex) > VSMALL)
            {
                oIndex = nIndex;

                scalar percent =
                (
                    100.0 -
                    (
                        (100.0 * mesh.stack(tIndex).size())
                      / (stackSize + VSMALL)
                    )
                );

                Info << "  Swap Progress: " << percent << "% :"
                     << "  Total: " << mesh.status(1)
                     << "             \r"
                     << flush;

                reported = true;
            }
        }

        // Retrieve the index for this face
        label fIndex = mesh.stack(tIndex).pop();

        // Perform a Delaunay test and check if a flip is necesary.
        bool failed = mesh.testDelaunay(fIndex);

        if (failed)
        {
            if (thread->master())
            {
                // Swap this face.
                mesh.swapQuadFace(fIndex);
            }
            else
            {
                // Push this on to the master stack
                mesh.stack(0).push(fIndex);
            }
        }
    }

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::STOP);
    }

    if (reported)
    {
        Info << "  Swap Progress: 100% :"
             << "  Total: " << mesh.status(1)
             << "             \r"
             << endl;
    }
}


// 3D Edge-swapping engine
void dynamicTopoFvMesh::swap3DEdges
(
    void *argument
)
{
    // Recast the argument
    meshHandler *thread = static_cast<meshHandler*>(argument);

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::START);
    }

    dynamicTopoFvMesh& mesh = thread->reference();

    // Figure out which thread this is...
    label tIndex = mesh.self();

    // Dynamic programming variables
    labelList m;
    PtrList<scalarListList> Q;
    PtrList<labelListList> K, triangulations;

    // Allocate dynamic programming tables
    mesh.initTables(m, Q, K, triangulations);

    // Set the timer
    clockTime sTimer;

    bool reported = false;
    label stackSize = mesh.stack(tIndex).size();
    scalar interval = mesh.reportInterval(), oIndex = 0.0, nIndex = 0.0;

    oIndex = ::floor(sTimer.elapsedTime() / interval);

    // Pick edges off the stack
    while (!mesh.stack(tIndex).empty())
    {
        // Report progress
        if (thread->master())
        {
            // Update the index, if its changed
            nIndex = ::floor(sTimer.elapsedTime() / interval);

            if ((nIndex - oIndex) > VSMALL)
            {
                oIndex = nIndex;

                scalar percent =
                (
                    100.0 -
                    (
                        (100.0 * mesh.stack(tIndex).size())
                      / (stackSize + VSMALL)
                    )
                );

                Info << "  Swap Progress: " << percent << "% :"
                     << "  Surface: " << mesh.status(2)
                     << ", Total: " << mesh.status(1)
                     << "             \r"
                     << flush;

                reported = true;
            }
        }

        // Retrieve an edge from the stack
        label eIndex = mesh.stack(tIndex).pop();

        // Compute the minimum quality of cells around this edge
        scalar minQuality = mesh.computeMinQuality(eIndex);

        // Check if this edge is on a bounding curve
        if (mesh.checkBoundingCurve(eIndex))
        {
            continue;
        }

        // Fill the dynamic programming tables
        if (mesh.fillTables(eIndex, minQuality, m, Q, K, triangulations))
        {
            // Check if edge-swapping is required.
            if (mesh.checkQuality(eIndex, m, Q, minQuality))
            {
                if (thread->master())
                {
                    // Remove this edge according to the swap sequence
                    mesh.removeEdgeFlips(eIndex, minQuality, K, triangulations);
                }
                else
                {
                    // Push this on to the master stack
                    mesh.stack(0).push(eIndex);
                }
            }
        }
    }

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::STOP);
    }

    if (reported)
    {
        Info << "  Swap Progress: 100% :"
             << "  Surface: " << mesh.status(2)
             << ", Total: " << mesh.status(1)
             << "             \r"
             << endl;
    }
}


// Edge refinement engine
void dynamicTopoFvMesh::edgeRefinementEngine
(
    void *argument
)
{
    // Loop through all edges and bisect/collapse by the criterion:
    // Bisect when edge-length > ratioMax_*lengthScale
    // Collapse when edge-length < ratioMin_*lengthScale

    // Recast the argument
    meshHandler *thread = static_cast<meshHandler*>(argument);

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::START);
    }

    dynamicTopoFvMesh& mesh = thread->reference();

    // Figure out which thread this is...
    label tIndex = mesh.self();

    // Set the timer
    clockTime sTimer;

    bool reported = false;
    label stackSize = mesh.stack(tIndex).size();
    scalar interval = mesh.reportInterval(), oIndex = 0.0, nIndex = 0.0;

    oIndex = ::floor(sTimer.elapsedTime() / interval);

    while (!mesh.stack(tIndex).empty())
    {
        // Update the index, if its changed
        // Report progress
        if (thread->master())
        {
            nIndex = ::floor(sTimer.elapsedTime() / interval);

            if ((nIndex - oIndex) > VSMALL)
            {
                oIndex = nIndex;

                scalar percent =
                (
                    100.0 -
                    (
                        (100.0 * mesh.stack(tIndex).size())
                      / (stackSize + VSMALL)
                    )
                );

                Info << "  Refinement Progress: " << percent << "% :"
                     << "  Bisections: " << mesh.status(3)
                     << ", Collapses: " << mesh.status(4)
                     << ", Total: " << mesh.status(0)
                     << "             \r"
                     << flush;

                reported = true;
            }
        }

        // Retrieve an entity from the stack
        label eIndex = mesh.stack(tIndex).pop();

        if (mesh.checkBisection(eIndex))
        {
            if (thread->master())
            {
                // Bisect this edge
                mesh.bisectEdge(eIndex);
            }
            else
            {
                // Push this on to the master stack
                mesh.stack(0).push(eIndex);
            }
        }
        else
        if (mesh.checkCollapse(eIndex))
        {
            if (thread->master())
            {
                // Collapse this edge
                mesh.collapseEdge(eIndex);
            }
            else
            {
                // Push this on to the master stack
                mesh.stack(0).push(eIndex);
            }
        }
    }

    if (thread->slave())
    {
        thread->sendSignal(meshHandler::STOP);
    }

    if (reported)
    {
        Info << "  Refinement Progress: 100% :"
             << "  Bisections: " << mesh.status(3)
             << ", Collapses: " << mesh.status(4)
             << ", Total: " << mesh.status(0)
             << "             \r"
             << endl;
    }
}


// Remove 2D sliver cells from the mesh
void dynamicTopoFvMesh::remove2DSlivers()
{
    // Sort by sliver-quality.
    labelList cIndices(thresholdSlivers_.toc());
    SortableList<scalar> values(cIndices.size());

    // Fill-in values to sort by...
    forAll(cIndices, indexI)
    {
        values[indexI] = thresholdSlivers_[cIndices[indexI]];
    }

    // Explicitly sort by quality value.
    values.sort();

    const labelList& indices = values.indices();

    if (debug)
    {
        if (thresholdSlivers_.size())
        {
            Info << "Sliver list: " << endl;

            forAll(indices, indexI)
            {
                label cIndex = cIndices[indices[indexI]];

                Info << " Cell: " << cIndex
                     << " Quality: " << thresholdSlivers_[cIndex]
                     << endl;
            }

            if (debug > 1)
            {
                writeVTK("sliverCells", cIndices, 3);
            }
        }
    }

    forAll(indices, indexI)
    {
        // Fetch the cell
        label cIndex = cIndices[indices[indexI]];

        // Ensure that this cell actually exists.
        if (cells_[cIndex].empty())
        {
            continue;
        }

        const cell& cellToCheck = cells_[cIndex];

        // Find an appropriate quad-face
        label fIndex = -1;

        forAll(cellToCheck, faceI)
        {
            if (faces_[cellToCheck[faceI]].size() == 4)
            {
                fIndex = cellToCheck[faceI];
                break;
            }
        }

        // Find the prism faces
        FixedList<label,2> cBdyIndex(-1), cIntIndex(-1);
        FixedList<face,2> cBdyFace, cIntFace;

        meshOps::findPrismFaces
        (
            fIndex,
            cIndex,
            faces_,
            cells_,
            neighbour_,
            cBdyFace,
            cBdyIndex,
            cIntFace,
            cIntIndex
        );

        if (debug > 1)
        {
            InfoIn("void dynamicTopoFvMesh::remove2DSlivers()")
                << nl
                << " Considering cell: " << cIndex
                << ":: " << cells_[cIndex]
                << " for sliver removal."
                << " Using face: " << fIndex
                << ":: " << faces_[fIndex]
                << endl;
        }

        label triFace = cBdyIndex[0], triEdge = -1;

        // Find the common-edge between quad-tri faces
        meshOps::findCommonEdge
        (
            fIndex,
            triFace,
            faceEdges_,
            triEdge
        );

        // Find the isolated point.
        label ptIndex = -1, nextPtIndex = -1;

        const edge& edgeToCheck = edges_[triEdge];

        meshOps::findIsolatedPoint
        (
            faces_[triFace],
            edgeToCheck,
            ptIndex,
            nextPtIndex
        );

        // Determine the interior faces connected to each edge-point.
        label firstFace = -1, secondFace = -1;

        if (cIntFace[0].which(edgeToCheck[0]) > -1)
        {
            firstFace  = cIntIndex[0];
            secondFace = cIntIndex[1];
        }
        else
        {
            firstFace  = cIntIndex[1];
            secondFace = cIntIndex[0];
        }

        point ec =
        (
            linePointRef
            (
                points_[edgeToCheck[0]],
                points_[edgeToCheck[1]]
            ).centre()
        );

        FixedList<vector, 2> p(vector::zero), q(vector::zero);
        FixedList<scalar, 2> proj(0.0);

        // Find the projection on the edge.
        forAll(edgeToCheck, pointI)
        {
            p[pointI] = (points_[ptIndex] - points_[edgeToCheck[pointI]]);
            q[pointI] = (ec - points_[edgeToCheck[pointI]]);

            q[pointI] /= (mag(q[pointI]) + VSMALL);

            proj[pointI] = (p[pointI] & q[pointI]);
        }

        // Take action based on the magnitude of the projection.
        if (mag(proj[0]) < VSMALL)
        {
            if (collapseQuadFace(firstFace).type() > 0)
            {
                statistics_[7]++;
            }

            continue;
        }

        if (mag(proj[1]) < VSMALL)
        {
            if (collapseQuadFace(secondFace).type() > 0)
            {
                statistics_[7]++;
            }

            continue;
        }

        if (proj[0] > 0.0 && proj[1] < 0.0)
        {
            changeMap map = bisectQuadFace(firstFace, changeMap());

            // Loop through added faces, and collapse
            // the appropriate one
            const List<FixedList<label,2> >& aF = map.addedFaceList();

            forAll(aF, faceI)
            {
                if
                (
                    (owner_[aF[faceI][0]] == cIndex) &&
                    (aF[faceI][0] != firstFace)
                )
                {
                    if (collapseQuadFace(aF[faceI][0]).type() > 0)
                    {
                        statistics_[7]++;
                    }

                    break;
                }
            }

            continue;
        }

        if (proj[0] < 0.0 && proj[1] > 0.0)
        {
            changeMap map = bisectQuadFace(secondFace, changeMap());

            // Loop through added faces, and collapse
            // the appropriate one
            const List<FixedList<label,2> >& aF = map.addedFaceList();

            forAll(aF, faceI)
            {
                if
                (
                    (owner_[aF[faceI][0]] == cIndex) &&
                    (aF[faceI][0] != secondFace)
                )
                {
                    if (collapseQuadFace(aF[faceI][0]).type() > 0)
                    {
                        statistics_[7]++;
                    }

                    break;
                }
            }

            continue;
        }

        if (proj[0] > 0.0 && proj[1] > 0.0)
        {
            changeMap map = bisectQuadFace(fIndex, changeMap());

            // Loop through added faces, and collapse
            // the appropriate one
            const List<FixedList<label,2> >& aF = map.addedFaceList();

            forAll(aF, faceI)
            {
                if
                (
                    (owner_[aF[faceI][0]] == cIndex) &&
                    (aF[faceI][0] != fIndex)
                )
                {
                    if (collapseQuadFace(aF[faceI][0]).type() > 0)
                    {
                        statistics_[7]++;
                    }

                    break;
                }
            }

            continue;
        }
    }

    // Clear out the list
    thresholdSlivers_.clear();
}


// Identify the sliver type in 3D
const changeMap dynamicTopoFvMesh::identifySliverType
(
    const label cIndex
) const
{
    changeMap map;

    // Ensure that this cell actually exists.
    if (cells_[cIndex].empty())
    {
        return map;
    }

    label fourthPoint = -1;
    scalar minDistance = GREAT;
    face tFace(3), testFace(3), faceToCheck(3);
    FixedList<edge, 2> edgeToCheck(edge(-1,-1));

    const cell& cellToCheck = cells_[cIndex];

    // Find the point-face pair with minimum perpendicular distance
    forAll(cellToCheck, faceI)
    {
        label isolatedPoint = -1;
        label nextFaceI = cellToCheck.fcIndex(faceI);

        // Pick two faces from this cell.
        const face& currFace = faces_[cellToCheck[faceI]];
        const face& nextFace = faces_[cellToCheck[nextFaceI]];

        // Get the fourth point
        forAll(nextFace, pointI)
        {
            if
            (
                nextFace[pointI] != currFace[0]
             && nextFace[pointI] != currFace[1]
             && nextFace[pointI] != currFace[2]
            )
            {
                isolatedPoint = nextFace[pointI];

                // Configure a triangular face with correct orientation.
                if (owner_[cellToCheck[faceI]] == cIndex)
                {
                    testFace[0] = currFace[2];
                    testFace[1] = currFace[1];
                    testFace[2] = currFace[0];
                }
                else
                {
                    testFace[0] = currFace[0];
                    testFace[1] = currFace[1];
                    testFace[2] = currFace[2];
                }

                break;
            }
        }

        // Obtain the unit normal.
        vector testNormal = testFace.normal(points_);

        testNormal /= (mag(testNormal) + VSMALL);

        // Project the isolated point onto the face.
        vector p = points_[isolatedPoint] - points_[testFace[0]];
        vector q = p - ((p & testNormal)*testNormal);

        // Compute the distance
        scalar distance = mag(p - q);

        // Is it the least so far?
        if (distance < minDistance)
        {
            // Use this point-face pair.
            fourthPoint = isolatedPoint;
            tFace = testFace;
            minDistance = distance;
        }
    }

    // Obtain the face-normal.
    vector refArea = tFace.normal(points_);

    // Normalize it.
    vector n = refArea/mag(refArea);

    // Define edge-vectors.
    vector r1 = points_[tFace[1]] - points_[tFace[0]];
    vector r2 = points_[tFace[2]] - points_[tFace[1]];
    vector r3 = points_[tFace[0]] - points_[tFace[2]];

    // Project the fourth point onto the face.
    vector r4 = points_[fourthPoint] - points_[tFace[0]];

    r4 = r4 - ((r4 & n)*n);

    // Define the two other vectors.
    vector r5 = r4 - r1;
    vector r6 = r5 - r2;

    // Calculate three signed triangle areas, using tFace[0] as the origin.
    scalar t1 = n & (0.5 * (r1 ^ r4));
    scalar t2 = n & (0.5 * (r2 ^ r5));
    scalar t3 = n & (0.5 * (r3 ^ r6));

    // Determine sliver types based on are magnitudes.
    if (t1 > 0 && t2 > 0 && t3 > 0)
    {
        // Region R0: Cap cell.
        map.type() = 2;
        map.apexPoint() = fourthPoint;

        faceToCheck[0] = tFace[0];
        faceToCheck[1] = tFace[1];
        faceToCheck[2] = tFace[2];
    }

    if (t1 < 0 && t2 > 0 && t3 > 0)
    {
        // Region R1: Sliver cell.
        map.type() = 1;

        edgeToCheck[0][0] = tFace[2];
        edgeToCheck[0][1] = fourthPoint;
        edgeToCheck[1][0] = tFace[0];
        edgeToCheck[1][1] = tFace[1];
    }

    if (t1 > 0 && t2 < 0 && t3 > 0)
    {
        // Region R2: Sliver cell.
        map.type() = 1;

        edgeToCheck[0][0] = tFace[0];
        edgeToCheck[0][1] = fourthPoint;
        edgeToCheck[1][0] = tFace[1];
        edgeToCheck[1][1] = tFace[2];
    }

    if (t1 > 0 && t2 > 0 && t3 < 0)
    {
        // Region R3: Sliver cell.
        map.type() = 1;

        edgeToCheck[0][0] = tFace[1];
        edgeToCheck[0][1] = fourthPoint;
        edgeToCheck[1][0] = tFace[2];
        edgeToCheck[1][1] = tFace[0];
    }

    if (t1 < 0 && t2 > 0 && t3 < 0)
    {
        // Region R4: Cap cell.
        map.type() = 2;
        map.apexPoint() = tFace[0];

        faceToCheck[0] = tFace[1];
        faceToCheck[1] = tFace[2];
        faceToCheck[2] = fourthPoint;
    }

    if (t1 < 0 && t2 < 0 && t3 > 0)
    {
        // Region R5: Cap cell.
        map.type() = 2;
        map.apexPoint() = tFace[1];

        faceToCheck[0] = tFace[2];
        faceToCheck[1] = tFace[0];
        faceToCheck[2] = fourthPoint;
    }

    if (t1 > 0 && t2 < 0 && t3 < 0)
    {
        // Region R6: Cap cell.
        map.type() = 2;
        map.apexPoint() = tFace[2];

        faceToCheck[0] = tFace[0];
        faceToCheck[1] = tFace[1];
        faceToCheck[2] = fourthPoint;
    }

    // See if an over-ride to wedge/spade is necessary.
    // Obtain a reference area magnitude.
    scalar refMag = 0.1*(refArea & n);

    if (mag(t1) < refMag)
    {
        if (mag(t3) < refMag)
        {
            // Wedge case: Too close to point [0]
            map.type() = 4;

            edgeToCheck[0][0] = tFace[0];
            edgeToCheck[0][1] = fourthPoint;
        }
        else
        if (mag(t2) < refMag)
        {
            // Wedge case: Too close to point [1]
            map.type() = 4;

            edgeToCheck[0][0] = tFace[1];
            edgeToCheck[0][1] = fourthPoint;
        }
        else
        if ((mag(t2) > refMag) && (mag(t3) > refMag))
        {
            // Spade case: Too close to edge vector r1
            map.type() = 3;
            map.apexPoint() = fourthPoint;

            edgeToCheck[0][0] = tFace[0];
            edgeToCheck[0][1] = tFace[1];
        }
    }

    if (mag(t2) < refMag)
    {
        if (mag(t3) < refMag)
        {
            // Wedge case: Too close to point [2]
            map.type() = 4;

            edgeToCheck[0][0] = tFace[2];
            edgeToCheck[0][1] = fourthPoint;
        }
        else
        if ((mag(t1) > refMag) && (mag(t3) > refMag))
        {
            // Spade case: Too close to edge vector r2
            map.type() = 3;
            map.apexPoint() = fourthPoint;

            edgeToCheck[0][0] = tFace[1];
            edgeToCheck[0][1] = tFace[2];
        }
    }

    if (mag(t3) < refMag)
    {
        if ((mag(t1) > refMag) && (mag(t2) > refMag))
        {
            // Spade case: Too close to edge vector r3
            map.type() = 3;
            map.apexPoint() = fourthPoint;

            edgeToCheck[0][0] = tFace[2];
            edgeToCheck[0][1] = tFace[0];
        }
    }

    // Determine appropriate information for sliver exudation.
    switch (map.type())
    {
        case 1:
        {
            FixedList<bool, 2> foundEdge(false);

            // Search the cell-faces for first and second edges.
            forAll(cellToCheck, faceI)
            {
                const labelList& fEdges = faceEdges_[cellToCheck[faceI]];

                forAll(fEdges, edgeI)
                {
                    const edge& thisEdge = edges_[fEdges[edgeI]];

                    if (thisEdge == edgeToCheck[0])
                    {
                        map.firstEdge() = fEdges[edgeI];

                        foundEdge[0] = true;
                    }

                    if (thisEdge == edgeToCheck[1])
                    {
                        map.secondEdge() = fEdges[edgeI];

                        foundEdge[1] = true;
                    }
                }

                if (foundEdge[0] && foundEdge[1])
                {
                    break;
                }
            }

            break;
        }

        case 2:
        {
            // Search the cell-faces for opposing faces.
            forAll(cellToCheck, faceI)
            {
                const face& thisFace = faces_[cellToCheck[faceI]];

                if (triFace::compare(triFace(thisFace), triFace(faceToCheck)))
                {
                    map.opposingFace() = cellToCheck[faceI];

                    break;
                }
            }

            break;
        }

        case 3:
        case 4:
        {
            bool foundEdge = false;

            // Search the cell-faces for first edge.
            forAll(cellToCheck, faceI)
            {
                const labelList& fEdges = faceEdges_[cellToCheck[faceI]];

                forAll(fEdges, edgeI)
                {
                    const edge& thisEdge = edges_[fEdges[edgeI]];

                    if (thisEdge == edgeToCheck[0])
                    {
                        map.firstEdge() = fEdges[edgeI];

                        foundEdge = true;
                    }
                }

                if (foundEdge)
                {
                    break;
                }
            }

            break;
        }

        default:
        {
            WarningIn
            (
                "void dynamicTopoFvMesh::identifySliverType"
                "(const label cIndex) const"
            )
                << nl << "Could not identify sliver type for cell: "
                << cIndex
                << endl;
        }
    }

    if (debug > 2)
    {
        Pout << "Cell: " << cIndex
             << " Identified sliver type as: "
             << map.type() << endl;
    }

    // Return the result.
    return map;
}


// Remove sliver cells
void dynamicTopoFvMesh::removeSlivers()
{
    if (!edgeRefinement_)
    {
        return;
    }

    // Check if a removeSlivers entry was found in the dictionary
    if (dict_.subDict("dynamicTopoFvMesh").found("removeSlivers"))
    {
        Switch rs =
        (
            dict_.subDict("dynamicTopoFvMesh").lookup("removeSlivers")
        );

        if (!rs)
        {
            return;
        }
    }

    // Invoke the 2D sliver removal routine
    if (twoDMesh_)
    {
        remove2DSlivers();
        return;
    }

    // Sort by sliver-quality.
    labelList cIndices(thresholdSlivers_.toc());
    SortableList<scalar> values(cIndices.size());

    // Fill-in values to sort by...
    forAll(cIndices, indexI)
    {
        values[indexI] = thresholdSlivers_[cIndices[indexI]];
    }

    // Explicitly sort by quality value.
    values.sort();

    const labelList& indices = values.indices();

    if (debug && thresholdSlivers_.size())
    {
        Info << "Sliver list: " << endl;

        forAll(indices, indexI)
        {
            label cIndex = cIndices[indices[indexI]];

            Info << " Cell: " << cIndex
                 << " Quality: " << thresholdSlivers_[cIndex]
                 << endl;
        }

        if (debug > 1)
        {
            writeVTK("sliverCells", cIndices, 3);
        }
    }

    forAll(indices, indexI)
    {
        // Fetch the cell index
        label cIndex = cIndices[indices[indexI]];

        // Identify the sliver type.
        changeMap map = identifySliverType(cIndex);

        if (debug)
        {
            InfoIn("void dynamicTopoFvMesh::removeSlivers()")
                << nl << "Removing Cell: " << cIndex
                << " of sliver type: " << map.type()
                << " with quality: " << thresholdSlivers_[cIndex]
                << endl;
        }

        bool success = false;

        // Take action based on the type of sliver.
        switch (map.type())
        {
            case -1:
            {
                // Sliver cell was removed by a prior operation.
                // Nothing needs to be done.
                break;
            }

            case 1:
            {
                // Sliver cell.
                // Determine which edges need to be bisected.
                label firstEdge = map.firstEdge();
                label secondEdge = map.secondEdge();

                // Force bisection on both edges.
                changeMap firstMap  = bisectEdge(firstEdge, false, true);
                changeMap secondMap = bisectEdge(secondEdge, false, true);

                // Collapse the intermediate edge.
                // Since we don't know which edge it is, search
                // through recently added edges and compare.
                edge edgeToCheck
                (
                    firstMap.addedPointList()[0][0],
                    secondMap.addedPointList()[0][0]
                );

                bool foundCollapseEdge = false;

                const List<FixedList<label,2> >& firstMapEdges =
                (
                    firstMap.addedEdgeList()
                );

                const List<FixedList<label,2> >& secondMapEdges =
                (
                    secondMap.addedEdgeList()
                );

                // Loop through the first list.
                forAll(firstMapEdges, edgeI)
                {
                    const edge& thisEdge =
                    (
                        edges_[firstMapEdges[edgeI][0]]
                    );

                    if (thisEdge == edgeToCheck)
                    {
                        // Collapse this edge.
                        if
                        (
                            collapseEdge
                            (
                                firstMapEdges[edgeI][0],
                                -1,
                                false,
                                true
                            ).type() > 0
                        )
                        {
                            success = true;
                        }

                        foundCollapseEdge = true;
                        break;
                    }
                }

                // Loop through the second list.
                if (!foundCollapseEdge)
                {
                    forAll(secondMapEdges, edgeI)
                    {
                        const edge& thisEdge =
                        (
                            edges_[secondMapEdges[edgeI][0]]
                        );

                        if (thisEdge == edgeToCheck)
                        {
                            // Collapse this edge.
                            if
                            (
                                collapseEdge
                                (
                                    secondMapEdges[edgeI][0],
                                    -1,
                                    false,
                                    true
                                ).type() > 0
                            )
                            {
                                success = true;
                            }

                            break;
                        }
                    }
                }

                break;
            }

            case 2:
            {
                // Cap cell.
                label opposingFace = map.opposingFace();

                // Force trisection of the opposing face.
                changeMap faceMap =
                (
                    trisectFace(opposingFace, false, true)
                );

                // Collapse the intermediate edge.
                // Since we don't know which edge it is, search
                // through recently added edges and compare.
                edge edgeToCheck
                (
                    map.apexPoint(),
                    faceMap.addedPointList()[0][0]
                );

                const List<FixedList<label,2> >& faceMapEdges =
                (
                    faceMap.addedEdgeList()
                );

                forAll(faceMapEdges, edgeI)
                {
                    const edge& thisEdge = edges_[faceMapEdges[edgeI][0]];

                    if (thisEdge == edgeToCheck)
                    {
                        // Collapse this edge.
                        if
                        (
                            collapseEdge
                            (
                                faceMapEdges[edgeI][0],
                                -1,
                                false,
                                true
                            ).type() > 0
                        )
                        {
                            success = true;
                        }

                        break;
                    }
                }

                break;
            }

            case 3:
            {
                // Spade cell.

                // Force bisection on the first edge.
                changeMap firstMap =
                (
                    bisectEdge(map.firstEdge(), false, true)
                );

                // Collapse the intermediate edge.
                // Since we don't know which edge it is, search
                // through recently added edges and compare.
                edge edgeToCheck
                (
                    map.apexPoint(),
                    firstMap.addedPointList()[0][0]
                );

                const List<FixedList<label,2> >& firstMapEdges =
                (
                    firstMap.addedEdgeList()
                );

                // Loop through the first list.
                forAll(firstMapEdges, edgeI)
                {
                    const edge& thisEdge = edges_[firstMapEdges[edgeI][0]];

                    if (thisEdge == edgeToCheck)
                    {
                        // Collapse this edge.
                        if
                        (
                            collapseEdge
                            (
                                firstMapEdges[edgeI][0],
                                -1,
                                false,
                                true
                            ).type() > 0
                        )
                        {
                            success = true;
                        }

                        break;
                    }
                }

                break;
            }

            case 4:
            {
                // Wedge cell.

                // Collapse the first edge.
                if
                (
                    collapseEdge
                    (
                        map.firstEdge(),
                        -1,
                        false,
                        true
                    ).type() > 0
                )
                {
                    success = true;
                }

                break;
            }

            default:
            {
                WarningIn("void dynamicTopoFvMesh::removeSlivers()")
                    << nl << "Could not identify sliver type for cell: "
                    << cIndex
                    << endl;
            }
        }

        // Increment the count for successful sliver removal
        if (success)
        {
            statistics_[7]++;
        }
    }

    // Clear out the list
    thresholdSlivers_.clear();
}


// MultiThreaded topology modifier
void dynamicTopoFvMesh::threadedTopoModifier()
{
    // Remove sliver cells first.
    removeSlivers();

    // Coupled entities to avoid during normal modification
    labelHashSet entities;

    // Set the thread scheduling sequence
    labelList topoSequence(threader_->getNumThreads());

    // Linear sequence from 1 to nThreads
    forAll(topoSequence, indexI)
    {
        topoSequence[indexI] = indexI + 1;
    }

    if (edgeRefinement_)
    {
        // Initialize stacks
        initStacks(entities);

        // Execute threads
        if (threader_->multiThreaded())
        {
            executeThreads(topoSequence, handlerPtr_, &edgeRefinementEngine);
        }

        // Set the master thread to implement modifications
        edgeRefinementEngine(&(handlerPtr_[0]));

        if (debug)
        {
            Info << nl << "Edge Bisection/Collapse complete." << endl;
        }
    }

    // Re-Initialize stacks
    initStacks(entities);

    // Execute threads
    if (threader_->multiThreaded())
    {
        if (twoDMesh_)
        {
            executeThreads(topoSequence, handlerPtr_, &swap2DEdges);
        }
        else
        {
            executeThreads(topoSequence, handlerPtr_, &swap3DEdges);
        }
    }

    // Set the master thread to implement modifications
    if (twoDMesh_)
    {
        swap2DEdges(&(handlerPtr_[0]));
    }
    else
    {
        swap3DEdges(&(handlerPtr_[0]));
    }

    if (debug)
    {
        Info << nl << "Edge Swapping complete." << endl;
    }
}


// Reset the mesh and generate mapping information
//  - Return true if topology changes were made.
//  - Return false otherwise (motion only)
bool dynamicTopoFvMesh::resetMesh()
{
    // Reduce across processors.
    reduce(topoChangeFlag_, orOp<bool>());

    if (topoChangeFlag_)
    {
        // Write out statistics
        Info << " Bisections :: Total: " << status(3)
             << ", Surface: " << status(5) << endl;
        Info << " Collapses  :: Total: " << status(4)
             << ", Surface: " << status(6) << endl;
        Info << " Swaps      :: Total: " << status(1)
             << ", Surface: " << status(2) << endl;

        if (status(7))
        {
            Info << " Slivers    :: " << status(7) << endl;
        }

        // Set sizes for mapping
        faceWeights_.setSize(facesFromFaces_.size(), scalarField(0));
        faceCentres_.setSize(facesFromFaces_.size(), vectorField(0));
        cellWeights_.setSize(cellsFromCells_.size(), scalarField(0));
        cellCentres_.setSize(cellsFromCells_.size(), vectorField(0));

        // Fetch the match tolerance for mapping
        scalar matchTol = Foam::debug::tolerances("meshOpsMatchTol", 1e-4);

        // Determine if mapping is to be skipped
        // Optionally skip mapping for remeshing-only / pre-processing
        const dictionary& meshSubDict = dict_.subDict("dynamicTopoFvMesh");

        bool skipMapping = false;

        if (meshSubDict.found("skipMapping") || mandatory_)
        {
            skipMapping = readBool(meshSubDict.lookup("skipMapping"));
        }

        clockTime mappingTimer;

        // Compute mapping weights for modified entities
        threadedMapping(matchTol, skipMapping);

        // Print out stats
        Info << " Mapping time: " << mappingTimer.elapsedTime()
             << " s" << endl;

        // Obtain references to zones, if any
        pointZoneMesh& pointZones = polyMesh::pointZones();
        faceZoneMesh& faceZones = polyMesh::faceZones();
        cellZoneMesh& cellZones = polyMesh::cellZones();

        // Allocate temporary lists for mesh-reset
        pointField points(nPoints_);
        pointField preMotionPoints(nPoints_);
        edgeList edges(nEdges_);
        faceList faces(nFaces_);
        labelList owner(nFaces_);
        labelList neighbour(nInternalFaces_);
        labelListList faceEdges(nFaces_);
        labelListList edgeFaces(nEdges_);
        labelList oldPatchStarts(oldPatchStarts_);
        labelList oldPatchNMeshPoints(oldPatchNMeshPoints_);
        labelListList pointZoneMap(pointZones.size());
        labelListList faceZonePointMap(faceZones.size());
        labelListList faceZoneFaceMap(faceZones.size());
        labelListList cellZoneMap(cellZones.size());

        // Obtain faceZone point maps before reordering
        List<Map<label> > oldFaceZonePointMaps(faceZones.size());

        forAll(faceZones, fzI)
        {
            oldFaceZonePointMaps[fzI] = faceZones[fzI]().meshPointMap();
        }

        clockTime reOrderingTimer;

        // Reorder the mesh and obtain current topological information
        reOrderMesh
        (
            points,
            preMotionPoints,
            edges,
            faces,
            owner,
            neighbour,
            faceEdges,
            edgeFaces,
            pointZoneMap,
            faceZoneFaceMap,
            cellZoneMap
        );

        // Print out stats
        Info << " Reordering time: " << reOrderingTimer.elapsedTime()
             << " s" << endl;

        // Obtain the patch-point maps before resetting the mesh
        List<Map<label> > oldPatchPointMaps(boundaryMesh().size());

        forAll(oldPatchPointMaps, patchI)
        {
            oldPatchPointMaps[patchI] = boundaryMesh()[patchI].meshPointMap();
        }

        topoMapper& fieldMapper = mapper_();

        // Set information for the mapping stage, prior to mesh reset
        fieldMapper.storeMeshInformation();

        // Set weighting information.
        // This takes over the weight data.
        fieldMapper.setFaceWeights
        (
            xferMove(faceWeights_),
            xferMove(faceCentres_)
        );

        fieldMapper.setCellWeights
        (
            xferMove(cellWeights_),
            xferMove(cellCentres_)
        );

        // Reset the mesh with pre-motion points
        polyMesh::resetPrimitives
        (
            xferCopy(preMotionPoints),
            xferMove(faces),
            xferMove(owner),
            xferMove(neighbour),
            patchSizes_,
            patchStarts_
        );

        // Check the dictionary to determine whether
        // edge-connectivity is to be stored on disk.
        // This usually benefits restart-time for large
        // cases, at the expense of disk-space.
        bool storePrimitives = false;

        if (meshSubDict.found("storeEdgePrimitives") || mandatory_)
        {
            storePrimitives =
            (
                readBool(meshSubDict.lookup("storeEdgePrimitives"))
            );
        }

        // Reset the edge mesh
        eMeshPtr_->resetPrimitives
        (
            edges,
            faceEdges,
            edgeFaces,
            edgePatchSizes_,
            edgePatchStarts_,
            true,
            (time().outputTime() && storePrimitives)
        );

        // Generate mapping for points on boundary patches
        labelListList patchPointMap(boundaryMesh().size());

        for (label i = 0; i < boundaryMesh().size(); i++)
        {
            // Obtain new patch mesh points after reset.
            const labelList& meshPointLabels = boundaryMesh()[i].meshPoints();

            patchNMeshPoints_[i] = meshPointLabels.size();

            patchPointMap[i].setSize(meshPointLabels.size(), -1);

            forAll(meshPointLabels, pointI)
            {
                label oldIndex = pointMap_[meshPointLabels[pointI]];

                // Check if the point existed before...
                if (oldIndex > -1)
                {
                    // Look for the old position on this patch.
                    Map<label>::const_iterator oIter =
                    (
                        oldPatchPointMaps[i].find(oldIndex)
                    );

                    // Add an entry if the point was found
                    if (oIter != oldPatchPointMaps[i].end())
                    {
                        patchPointMap[i][pointI] = oIter();
                    }
                }
            }
        }

        // Generate mapping for faceZone points
        forAll(faceZones, fzI)
        {
            // Obtain new face zone mesh points after reset.
            const labelList& meshPointLabels = faceZones[fzI]().meshPoints();

            faceZonePointMap[fzI].setSize(meshPointLabels.size(), -1);

            forAll(meshPointLabels, pointI)
            {
                label oldIndex = pointMap_[meshPointLabels[pointI]];

                // Check if the point existed before...
                if (oldIndex > -1)
                {
                    // Look for the old position on this patch.
                    Map<label>::const_iterator oIter =
                    (
                        oldFaceZonePointMaps[fzI].find(oldIndex)
                    );

                    // Add an entry if the point was found
                    if (oIter != oldFaceZonePointMaps[fzI].end())
                    {
                        faceZonePointMap[fzI][pointI] = oIter();
                    }
                }
            }
        }

        // Generate new mesh mapping information
        mapPolyMesh mpm
        (
            (*this),
            nOldPoints_,
            nOldFaces_,
            nOldCells_,
            pointMap_,
            pointsFromPoints_,
            faceMap_,
            facesFromPoints_,
            facesFromEdges_,
            facesFromFaces_,
            cellMap_,
            cellsFromPoints_,
            cellsFromEdges_,
            cellsFromFaces_,
            cellsFromCells_,
            reversePointMap_,
            reverseFaceMap_,
            reverseCellMap_,
            flipFaces_,
            patchPointMap,
            pointZoneMap,
            faceZonePointMap,
            faceZoneFaceMap,
            cellZoneMap,
            preMotionPoints,
            oldPatchStarts,
            oldPatchNMeshPoints,
            true
        );

        // Update the underlying mesh, and map all related fields
        updateMesh(mpm);

        // Perform a dummy movePoints to force V0 creation, if necessary
        movePoints(mpm.preMotionPoints());

        // Reset old-volumes
        resetMotion();
        setV0();

        // Correct volume fluxes on the old mesh
        fieldMapper.correctFluxes();

        // Clear mapper after use
        fieldMapper.clear();

        // Now move mesh to new points and
        // compute correct mesh-fluxes.
        movePoints(points);

        // Update the mesh-motion solver
        if (motionSolver_.valid())
        {
            motionSolver_->updateMesh(mpm);
        }

        // Clear unwanted member data
        addedFacePatches_.clear();
        addedEdgePatches_.clear();
        addedPointZones_.clear();
        addedFaceZones_.clear();
        addedCellZones_.clear();
        faceParents_.clear();
        cellParents_.clear();

        // Clear the deleted entity map
        deletedPoints_.clear();
        deletedEdges_.clear();
        deletedFaces_.clear();
        deletedCells_.clear();

        // Clear flipFaces / modPoints
        flipFaces_.clear();
        modPoints_.clear();

        // Set new sizes for the reverse maps
        reversePointMap_.setSize(nPoints_, -7);
        reverseEdgeMap_.setSize(nEdges_, -7);
        reverseFaceMap_.setSize(nFaces_, -7);
        reverseCellMap_.setSize(nCells_, -7);

        // Update "old" information
        nOldPoints_ = nPoints_;
        nOldEdges_ = nEdges_;
        nOldFaces_ = nFaces_;
        nOldCells_ = nCells_;
        nOldInternalFaces_ = nInternalFaces_;
        nOldInternalEdges_ = nInternalEdges_;

        for (label i = 0; i < boundaryMesh().size(); i++)
        {
            oldPatchSizes_[i] = patchSizes_[i];
            oldPatchStarts_[i] = patchStarts_[i];
            oldEdgePatchSizes_[i] = edgePatchSizes_[i];
            oldEdgePatchStarts_[i] = edgePatchStarts_[i];
            oldPatchNMeshPoints_[i] = patchNMeshPoints_[i];
        }

        // Basic checks for mesh-validity
        if (debug > 2)
        {
            checkMesh(true);
        }

        // Reset statistics
        statistics_ = 0;
    }
    else
    {
        // No topology changes were made.
        // Only execute mesh-motion.
        if (motionSolver_.valid())
        {
            movePoints(motionSolver_->curPoints());
        }
    }

    // Obtain mesh stats after topo-changes
    meshQuality(true);

    // Dump length-scale to disk, if requested.
    calculateLengthScale(true);

    // Reset and return flag
    if (topoChangeFlag_)
    {
        topoChangeFlag_ = false;
        return true;
    }

    // No changes were made.
    return false;
}


// Map all fields in time using a customized mapper
void dynamicTopoFvMesh::mapFields(const mapPolyMesh& meshMap) const
{
    if (debug)
    {
        Info << "void dynamicTopoFvMesh::mapFields(const mapPolyMesh&): "
             << "Mapping fv fields."
             << endl;
    }

    const topoMapper& fieldMapper = mapper_();

    // Set the mapPolyMesh object in the mapper
    fieldMapper.setMapper(meshMap);

    // Map all the volFields in the objectRegistry
    MapGeometricFields<scalar,fvPatchField,topoMapper,volMesh>
        (fieldMapper);
    MapGeometricFields<vector,fvPatchField,topoMapper,volMesh>
        (fieldMapper);
    MapGeometricFields<sphericalTensor,fvPatchField,topoMapper,volMesh>
        (fieldMapper);
    MapGeometricFields<symmTensor,fvPatchField,topoMapper,volMesh>
        (fieldMapper);
    MapGeometricFields<tensor,fvPatchField,topoMapper,volMesh>
        (fieldMapper);

    // Map all the surfaceFields in the objectRegistry
    MapGeometricFields<scalar,fvsPatchField,topoMapper,surfaceMesh>
        (fieldMapper);
    MapGeometricFields<vector,fvsPatchField,topoMapper,surfaceMesh>
        (fieldMapper);
    MapGeometricFields<sphericalTensor,fvsPatchField,topoMapper,surfaceMesh>
        (fieldMapper);
    MapGeometricFields<symmTensor,fvsPatchField,topoMapper,surfaceMesh>
        (fieldMapper);
    MapGeometricFields<tensor,fvsPatchField,topoMapper,surfaceMesh>
        (fieldMapper);
}


// Update the mesh for motion / topology changes.
//  - Return true if topology changes have occurred
bool dynamicTopoFvMesh::update()
{
    // Re-read options, in case they have been modified at run-time
    readOptionalParameters(true);

    // Set old point positions
    oldPoints_ = polyMesh::points();

    // Invoke mesh-motion solver and store new points
    if (motionSolver_.valid())
    {
        points_ = motionSolver_->newPoints()();
    }

    // Obtain mesh stats before topo-changes
    bool noSlivers = meshQuality(true);

    // Return if the interval is invalid,
    // not at re-mesh interval, or slivers are absent.
    // Handy while using only mesh-motion.
    if (interval_ < 0 || ((time().timeIndex() % interval_ != 0) && noSlivers))
    {
        return resetMesh();
    }

    // Calculate the edge length-scale for the mesh
    calculateLengthScale();

    // Track mesh topology modification time
    clockTime topoTimer;

    // Invoke the threaded topoModifier
    threadedTopoModifier();

    Info << " Topo modifier time: "
         << topoTimer.elapsedTime() << " s"
         << endl;

    // Apply all topology changes (if any) and reset mesh.
    return resetMesh();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void dynamicTopoFvMesh::operator=(const dynamicTopoFvMesh& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "dynamicTopoFvMesh::operator=(const dynamicTopoFvMesh&)"
        )
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


} // End namespace Foam

// ************************************************************************* //
