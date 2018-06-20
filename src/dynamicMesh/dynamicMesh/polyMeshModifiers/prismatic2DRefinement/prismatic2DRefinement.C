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

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

Notes
    Specialisation of polyhedralRefinement for 2D simulations on arbitrary
    prismatic meshes (tet, hex, etc).

\*---------------------------------------------------------------------------*/

#include "prismatic2DRefinement.H"
#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "foamTime.H"
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "meshTools.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "mapPolyMesh.H"
#include "emptyPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(prismatic2DRefinement, 0);
    addToRunTimeSelectionTable
    (
        polyMeshModifier,
        prismatic2DRefinement,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::prismatic2DRefinement::getAnchorLevel
(
    const label faceI,
    const label nPoints
) const
{
    // Get the face
    const face& f = mesh_.faces()[faceI];

    // Sanity check for expected number of points
    if (nPoints != 3 && nPoints != 4)
    {
        FatalErrorIn
        (
            "label prismatic2DRefinement::getAnchorLevel"
            "\n("
            "\n    const label faceI,"
            "\n    const label nPoints"
            "\n) const"
        )   << "Trying to find anchor level with " << nPoints << " points"
            << " smaller than anchor level. Only nPoints = 3 and 4 are"
            << " supported."
            << abort(FatalError);
    }

    // Sanity check: if we are expecting to find 4 points on a face which is not
    // on empty patch and we find a face with 3 points, issue an error
    if (f.size() <= 3 && nPoints == 4)
    {
        FatalErrorIn
        (
            "label prismatic2DRefinement::getAnchorLevel"
            "\n("
            "\n    const label faceI,"
            "\n    const label nPoints"
            "\n) const"
        )   << "Expected to find at least 4 points with level lower than"
            << " anchor level."
            << nl
            << "Make sure to call this function with nPoints = 4 only if"
            << " you do not expect triangular faces."
            << abort(FatalError);
    }

    if (f.size() <= nPoints)
    {
        return pointLevel_[f[findMaxLevel(f)]];
    }
    else
    {
        const label& ownLevel = cellLevel_[mesh_.faceOwner()[faceI]];

        if (countAnchors(f, ownLevel) >= nPoints)
        {
            return ownLevel;
        }
        else if (countAnchors(f, ownLevel + 1) >= nPoints)
        {
            return ownLevel + 1;
        }
        else
        {
            return -1;
        }
    }
}


void Foam::prismatic2DRefinement::calcLevel0EdgeLength()
{
    if (cellLevel_.size() != mesh_.nCells())
    {
        FatalErrorIn("scalar prismatic2DRefinement::getLevel0EdgeLength() const")
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size()
            << endl
            << "This might be because of a restart with inconsistent cellLevel."
            << abort(FatalError);
    }

    // Determine minimum edge length per refinement level
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar GREAT2 = sqr(GREAT);

    const label nLevels = gMax(cellLevel_) + 1;

    // Initialise typical edge length squared to dummy large value
    scalarField typEdgeLenSqr(nLevels, GREAT2);


    // 1. Look only at edges surrounded by cellLevel cells only
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Edge list containing cellLevel of connected cells.
    // Note:
    // - If not set = -1
    // - Levels of connected cells = from 0 to nLevels - 1
    // - If different levels = labelMax
    labelList edgeLevel(mesh_.nEdges(), -1);

    // Mark edges on empty patches. This will be used to skip edges on internal
    // faces and ordinary boundaries
    boolList edgeOnEmptyPatch(mesh_.nEdges() , false);

    // Get boundary mesh
    const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();

    // Get face-edge addressing
    const labelListList& meshFaceEdges = mesh_.faceEdges();

    // Loop through all patches
    forAll (boundaryMesh, patchI)
    {
        // Get current patch
        const polyPatch& curPatch = boundaryMesh[patchI];

        // Check whether this patch is emptyPolyPatch
        if (isA<emptyPolyPatch>(curPatch))
        {
            // Get start and end face labels
            const label startFaceI = curPatch.start();
            const label endFaceI = startFaceI + curPatch.size();

            // Mark all the faces and edges on the patch
            for (label faceI = startFaceI; faceI < endFaceI; ++faceI)
            {
                // Get edges of this face
                const labelList& curEdges = meshFaceEdges[faceI];

                // Mark all edges
                forAll (curEdges, i)
                {
                    edgeOnEmptyPatch[curEdges[i]] = true;
                }
            } // End for all patch faces
        } // End if the patch is empty
    } // End for all patches

    // Get cell edges
    const labelListList& meshCellEdges = mesh_.cellEdges();

    // Loop through cells
    forAll(cellLevel_, cellI)
    {
        // Get cell level and cell edges
        const label& cLevel = cellLevel_[cellI];
        const labelList& edges = meshCellEdges[cellI];

        // Loop through all edges
        forAll(edges, i)
        {
            const label& edgeI = edges[i];

            // Check whether the edge is on empty patch
            if (edgeOnEmptyPatch[edgeI])
            {
                if (edgeLevel[edgeI] == -1)
                {
                    // Edge level not set, mark it with current cell level
                    edgeLevel[edgeI] = cLevel;
                }
                else if (edgeLevel[edgeI] == labelMax)
                {
                    // Already marked as on different cellLevels
                }
                else if (edgeLevel[edgeI] != cLevel)
                {
                    // Edge level with different cell levels
                    edgeLevel[edgeI] = labelMax;
                }
            }
            // Else leave edgeLevel to dummy -1 value
        }
    }

    // Make sure that edges with different levels on different processors
    // are also marked. Do the same test (edgeLevel != cLevel) on coupled
    // edges.
    syncTools::syncEdgeList
    (
        mesh_,
        edgeLevel,
        ifEqEqAssignFirstOp<label, labelMax>(),
        labelMin,
        false               // no separation
    );

    // Now use the edgeLevel with a valid value to determine the
    // length per level.

    // Get mesh edges and points
    const edgeList& meshEdges = mesh_.edges();
    const pointField& meshPoints = mesh_.points();

    forAll(edgeLevel, edgeI)
    {
        // Get edge level
        const label& eLevel = edgeLevel[edgeI];

        // Edge has unique cell level (both cells have the same level)
        if (eLevel > -1 && eLevel < labelMax)
        {
            // Get the edge and calculate the square of the length
            const edge& e = meshEdges[edgeI];
            const scalar edgeLenSqr = magSqr(e.vec(meshPoints));

            typEdgeLenSqr[eLevel] = min(typEdgeLenSqr[eLevel], edgeLenSqr);
        }
        // Else typical edge length squared is not affected
    }

    // Get the minimum per level over all processors. Note: using minEqOp
    // because if the cells are not cubic, we end up using the smallest edge
    Pstream::listCombineGather(typEdgeLenSqr, minEqOp<scalar>());
    Pstream::listCombineScatter(typEdgeLenSqr);

    if (debug)
    {
        Pout<< "prismatic2DRefinement::calcLevel0EdgeLength() :"
            << " Initial edge lengths squared per refinementlevel:"
            << typEdgeLenSqr << endl;
    }


    // 2. For any levels where we haven't determined a valid length yet
    //    use any surrounding cell level. Here we use the max so we don't
    //    pick up levels between cellLevel and higher cellLevel (will have
    //    edges sized according to highest cellLevel)
    //    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialise maximum edge length to dummy large value
    scalarField maxEdgeLenSqr(nLevels, -GREAT2);

    // Loop through cells
    forAll(cellLevel_, cellI)
    {
        // Get cell level and cell edges
        const label& cLevel = cellLevel_[cellI];
        const labelList& cEdges = meshCellEdges[cellI];

        // Loop through edges
        forAll(cEdges, i)
        {
            // Get the edge
            const label& edgeI = cEdges[i];

            if (edgeOnEmptyPatch[edgeI])
            {
                // Edge is on empty patch, calculate the max edge length
                const edge& e = meshEdges[cEdges[i]];
                const scalar edgeLenSqr = magSqr(e.vec(meshPoints));

                maxEdgeLenSqr[cLevel] = max(maxEdgeLenSqr[cLevel], edgeLenSqr);
            }
            // Else typical max edge length squared is not affected
        }
    }

    // Get maximum per level over all processors
    Pstream::listCombineGather(maxEdgeLenSqr, maxEqOp<scalar>());
    Pstream::listCombineScatter(maxEdgeLenSqr);

    if (debug)
    {
        Pout<< "prismatic2DRefinement::calcLevel0EdgeLength() :"
            << " Basic edge lengths squared per refinementlevel:"
            << maxEdgeLenSqr << endl;
    }


    // 3. Combine the two sets of lengths
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(typEdgeLenSqr, levelI)
    {
        if
        (
            equal(typEdgeLenSqr[levelI], GREAT2)
         && maxEdgeLenSqr[levelI] >= 0
        )
        {
            typEdgeLenSqr[levelI] = maxEdgeLenSqr[levelI];
        }
    }

    if (debug)
    {
        Pout<< "prismatic2DRefinement::calcLevel0EdgeLength() :"
            << " Final edge lengths squared per refinementlevel:"
            << typEdgeLenSqr << endl;
    }

    // Find lowest length present across all levels
    level0EdgeLength_ = -1;

    forAll(typEdgeLenSqr, levelI)
    {
        const scalar& lenSqr = typEdgeLenSqr[levelI];

        if (lenSqr < GREAT2)
        {
            // Note: use power instead of left shift operator to multiply the
            // length with appropriate level. Easier to read. VV, 5/Jan/2018.
            level0EdgeLength_ = Foam::sqrt(lenSqr)*(Foam::pow(2, levelI));

            if (debug)
            {
                Pout<< "prismatic2DRefinement::calcLevel0EdgeLength() :"
                    << " Edge length: " << Foam::sqrt(lenSqr) << ", "
                    << "with edge level: " << levelI << ", "
                    << "has equivalent level0 lenght of:" << level0EdgeLength_
                    << endl;
            }

            break;
        }
    }

    if (level0EdgeLength_ == -1)
    {
        FatalErrorIn("prismatic2DRefinement::calcLevel0EdgeLength()")
            << "Problem in definition of typical edge length squared: "
            << typEdgeLenSqr << abort(FatalError);
    }
}


void Foam::prismatic2DRefinement::setInstance(const fileName& inst) const
{
    if (debug)
    {
        Pout<< "prismatic2DRefinement::setInstance(const fileName& inst)"
            << nl
            << "Resetting file instance of refinement data to " << inst
            << endl;
    }

    cellLevel_.instance() = inst;
    pointLevel_.instance() = inst;
}


void Foam::prismatic2DRefinement::extendMarkedCellsAcrossFaces
(
    boolList& markedCell
) const
{
    // Mark all faces for all marked cells
    const label nFaces = mesh_.nFaces();
    boolList markedFace(nFaces, false);

    // Get mesh cells
    const cellList& meshCells = mesh_.cells();

    // Loop through all cells
    forAll (markedCell, cellI)
    {
        if (markedCell[cellI])
        {
            // This cell is marked, get its faces
            const cell& cFaces = meshCells[cellI];

            forAll (cFaces, i)
            {
                markedFace[cFaces[i]] = true;
            }
        }
    }

    // Snyc the face list across processor boundaries
    syncTools::syncFaceList(mesh_, markedFace, orEqOp<bool>(), false);

    // Get necessary mesh data
    const label nInternalFaces = mesh_.nInternalFaces();
    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    // Internal faces
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        if (markedFace[faceI])
        {
            // Face will be marked, mark both owner and neighbour
            markedCell[owner[faceI]] = true;
            markedCell[neighbour[faceI]] = true;
        }
    }

    // Boundary faces
    for (label faceI = nInternalFaces; faceI < nFaces; ++faceI)
    {
        if (markedFace[faceI])
        {
            // Face will be markedd, mark owner
            markedCell[owner[faceI]] = true;
        }
    }
}


void Foam::prismatic2DRefinement::extendMarkedCellsAcrossPoints
(
    boolList& markedCell
) const
{
    // Mark all points for all marked cells
    const label nPoints = mesh_.nPoints();
    boolList markedPoint(nPoints, false);

    // Get cell points
    const labelListList& meshCellPoints = mesh_.cellPoints();

    // Loop through all cells
    forAll (markedCell, cellI)
    {
        if (markedCell[cellI])
        {
            // This cell is marked, get its points
            const labelList& cPoints = meshCellPoints[cellI];

            forAll (cPoints, i)
            {
                markedPoint[cPoints[i]] = true;
            }
        }
    }

    // Snyc point list across processor boundaries
    syncTools::syncPointList
    (
        mesh_,
        markedPoint,
        orEqOp<bool>(),
        true, // Default value
        true  // Apply separation for parallel cyclics
    );

    // Get point cells
    const labelListList& meshPointCells = mesh_.pointCells();

    // Loop through all points
    forAll (markedPoint, pointI)
    {
        if (markedPoint[pointI])
        {
            // This point is marked, mark all of its cells
            const labelList& pCells = meshPointCells[pointI];

            forAll (pCells, i)
            {
                markedCell[pCells[i]] = true;
            }
        }
    }
}


void Foam::prismatic2DRefinement::appendFaceSplitInfo
(
    const label& faceI,
    const boolList& edgeOnEmptyPatch,
    const labelList& edgeMidPoint,
    DynamicList<label>& splitFacesIntoTwo,
    DynamicList<Pair<label> >& splitFacesEmptyEdges
) const
{
    // First append the face into the list
    splitFacesIntoTwo.append(faceI);

    // Grab all edges of the face
    const labelList& curEdges = mesh_.faceEdges()[faceI];

    // Create placeholders for two edge levels. Initialise with -1
    // for sanity checks later on
    label emptyPatchEdgeI = -1;
    label emptyPatchEdgeJ = -1;

    // Count number of edges for sanity checks
    label nEdgesOnEmptyPatch = 0;

    // Collect the two edge labels found on the empty patch
    forAll (curEdges, i)
    {
        // Get edge index
        const label edgeI = curEdges[i];

        if (edgeOnEmptyPatch[edgeI])
        {
            // This edge is on empty patch, check whether its first,
            // second or invalid
            switch (nEdgesOnEmptyPatch)
            {
                case 0:
                    emptyPatchEdgeI = edgeI;
                    break;
                case 1:
                    emptyPatchEdgeJ = edgeI;
                    break;
                default:
                    FatalErrorIn
                    (
                        "prismatic2DRefinement::appendFaceSplitInfo(...)"
                    )   << "Found more than two edges on face " << faceI
                        << " on the empty patch."
                        << nl
                        << "Either this is not a valid 2D mesh or"
                        << " we are visiting wrong faces."
                        << abort(FatalError);
            }

            // Increment the counter
            ++nEdgesOnEmptyPatch;

        } // End if edge on empty patch

    } // End loop over all edges

    // Debug: additional check whether the two edges are marked for
    // refinement (they should be)
    if (debug)
    {
        if (edgeMidPoint[emptyPatchEdgeI] == -1)
        {
            FatalErrorIn
            (
                "prismatic2DRefinement::appendFaceSplitInfo(...)"
            )   << "Empty patch edge with index: " << emptyPatchEdgeI
                << " not marked for splitting"
                << nl
                << "Check edgeMidPoint selection algorithm."
                << abort(FatalError);
        }

        if (edgeMidPoint[emptyPatchEdgeI] == -1)
        {
            FatalErrorIn
            (
                "prismatic2DRefinement::appendFaceSplitInfo(...)"
            )   << "Empty patch edge with index: " << emptyPatchEdgeJ
                << " not marked for splitting"
                << nl
                << "Check edgeMidPoint selection algorithm."
                << abort(FatalError);
        }
    }

    // At this point, we should have the two edges we were looking
    // for, collect them into the list with additional sanity check
    if (emptyPatchEdgeI > -1 && emptyPatchEdgeJ > -1)
    {
        // Append the list of two edges in increasing order (just in
        // case we end up needing this information
        if (emptyPatchEdgeI < emptyPatchEdgeJ)
        {
            splitFacesEmptyEdges.append
            (
                Pair<label>(emptyPatchEdgeI, emptyPatchEdgeJ)
            );
        }
        else
        {
            splitFacesEmptyEdges.append
            (
                Pair<label>(emptyPatchEdgeJ, emptyPatchEdgeI)
            );
        }
    }
    else
    {
        FatalErrorIn
        (
            "prismatic2DRefinement::appendFaceSplitInfo(...)"
        )   << "Found invalid indices for edges on empty patches:"
            << nl
            << "emptyPatchEdgeI: " << emptyPatchEdgeI
            << ", emptyPatchEdgeJ: " << emptyPatchEdgeJ
            << nl
            << "Something went wrong. Check face edges."
            << abort(FatalError);
    }
}


void Foam::prismatic2DRefinement::setPrismatic2DRefinement
(
    polyTopoChange& ref
) const
{
    // Note: assumes that cellsToRefine_ are set prior to the function call

    // Reset refinementLevelIndicator field. Note: the list is cleared in
    // updateMesh member function after updating cell and point levels
    if (refinementLevelIndicator_.empty())
    {
        // List has been reset correctly, initialise it for this iteration
        refinementLevelIndicator_.setSize(mesh_.nCells(), UNCHANGED);
    }
    else
    {
        // List has not been reset correctly, issue an error
        FatalErrorIn
        (
            "prismatic2DRefinement::setPrismatic2DRefinement(...)"
        )   << "Refinement level indicator list has not been"
            << " reset properly." << nl
            << "Either the call to updateMesh() after performing"
            << " refinement has not been made or the call to"
            << " setPrismatic2DRefinement(...) and"
            << " setPrismatic2DUnrefinement(...) has not been made in"
            << " correct order." << nl
            << "Make sure to set refinement, set unrefinement and call"
            << " updateMesh after performing the topo change."
            << abort(FatalError);
    }

    if (debug)
    {
        // Write cells to refine into cell set
        const cellSet cSet
        (
            mesh_,
            "cellsToRefine",
            labelHashSet(cellsToRefine_)
        );

        cSet.write();

        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement(...)" << nl
            << "Writing " << cSet.size() << " cells to refine into cell set: "
            << cSet.objectPath()
            << endl;
    }


    // PART 1: Mark cells for refinement

    // Bool list that marks cells which will be refined
    boolList refineCellsMask(mesh_.nCells(), false);
    forAll (cellsToRefine_, i)
    {
        // Simply mark the cell as refined, there are no additional points to
        // add (in cell centre for example)
        refineCellsMask[cellsToRefine_[i]] = true;
    }

    // Write out cells to refine as a cell set for debug
    if (debug)
    {
        // Note: cellSet is actually a hash table of labels
        cellSet splitCells(mesh_, "splitCells", cellsToRefine_.size());

        forAll(refineCellsMask, cellI)
        {
            if (refineCellsMask[cellI])
            {
                // Cell is marked for refinement, insert into cellSet
                splitCells.insert(cellI);
            }
        }

        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement(...)" << nl
            << "Writing " << splitCells.size()
            << " cells to split to cellSet " << splitCells.objectPath()
            << endl;

        splitCells.write();
    }


    // PART 2: Mark edges for refinement and add points to edge centres

    if (debug)
    {
        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement(...)" << nl
            << "Allocating edge midpoints."
            << endl;
    }

    // First mark faces and edges on empty patches. This data is used in this
    // PART 2 (collecting edges) and also PART 3 (collecting faces)
    boolList faceOnEmptyPatch(mesh_.nFaces(), false);
    boolList edgeOnEmptyPatch(mesh_.nEdges(), false);

    // Get boundary
    const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();

    // Get face-edge addressing (for each face, a list of edges)
    const labelListList& meshFaceEdges = mesh_.faceEdges();

    // Loop through all patches
    forAll (boundaryMesh, patchI)
    {
        // Get current patch
        const polyPatch& curPatch = boundaryMesh[patchI];

        // Check whether this patch is emptyPolyPatch
        if (isA<emptyPolyPatch>(curPatch))
        {
            // Get start and end face labels
            const label startFaceI = curPatch.start();
            const label endFaceI = startFaceI + curPatch.size();

            // Mark all the faces and edges on the patch
            for (label faceI = startFaceI; faceI < endFaceI; ++faceI)
            {
                // Mark face
                faceOnEmptyPatch[faceI] = true;

                // Get edges of this face
                const labelList& curEdges = meshFaceEdges[faceI];

                // Mark all edges
                forAll (curEdges, i)
                {
                    edgeOnEmptyPatch[curEdges[i]] = true;
                }
            }
        }
    }

    // Now that we have marked faces and edges on empty patches, let's collect
    // refined edges. Refined edges are defined by having both their point
    // levels <= cell level, i.e. if any cell that gets split uses this edge
    // and the edge is on empty patch, the edge needs to be split

    // Get necessary mesh data
    const labelListList& meshCellEdges = mesh_.cellEdges();
    const edgeList& meshEdges = mesh_.edges();

    // Mid points for refined edge:
    // No need to split edge = -1
    // Label of introduced mid point > -1
    labelList edgeMidPoint(mesh_.nEdges(), -1);

    // Note: Loop over refined cells
    forAll (cellsToRefine_, i)
    {
        // Get cell index
        const label& cellI = cellsToRefine_[i];

        // Get edges of this cell
        const labelList& cEdges = meshCellEdges[cellI];

        forAll (cEdges, j)
        {
            // Get edge index and edge
            const label& edgeI = cEdges[j];
            const edge& e = meshEdges[edgeI];

            if (edgeOnEmptyPatch[edgeI])
            {
                // Edge is on empty patch, check whether it needs to be split
                if
                (
                    pointLevel_[e[0]] <= cellLevel_[cellI]
                 && pointLevel_[e[1]] <= cellLevel_[cellI]
                )
                {
                    // Point levels of both edge points are <= cell level, mark
                    // edge for splitting
                    edgeMidPoint[edgeI] = 12345;
                }
            }
            // Else nothing to do: can't split edges that are not on empty patch
        } // End for all edges of the refined cell
    } // End for all cells

    // Synchronize edgeMidPoint across coupled patches. Note: use max so that
    // any split takes precedence.
    syncTools::syncEdgeList
    (
        mesh_,
        edgeMidPoint,
        maxEqOp<label>(),
        labelMin,
        false               // no separation
    );

    // Now that the refinement trigger is synced, introduce edge points

    // Get necessary mesh data
    const pointField& meshPoints = mesh_.points();

    // Memory management
    {
        // Phase 1: calculate midpoints and sync. This is necessary if we don't
        // use binary format for writing and we slowly get differences.

        // Allocate storage for edge points
        pointField edgeMids(mesh_.nEdges(), point(-GREAT, -GREAT, -GREAT));

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] > -1)
            {
                // Edge marked to be split. Get edge centre
                edgeMids[edgeI] = meshEdges[edgeI].centre(meshPoints);
            }
        }

        // Sync across processor boundaries
        syncTools::syncEdgeList
        (
            mesh_,
            edgeMids,
            maxEqOp<vector>(),
            point(-GREAT, -GREAT, -GREAT),
            true                           // apply separation
        );

        // Phase 2: introduce points at the synced locations.
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] > -1)
            {
                edgeMidPoint[edgeI] = ref.setAction
                (
                    polyAddPoint
                    (
                        edgeMids[edgeI], // Point
                        -1,              // Appended point, no master ID
                        -1,              // Zone for point
                        true             // Supports a cell
                    )
                );
            }
        }
    } // End memory management for syncing/adding edge points

    // Write out edge mid points for split edges for debugging
    if (debug)
    {
        OFstream str(mesh_.time().path()/"edgeMidPoint.obj");

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] > -1)
            {
                // Get edge and write its cell centre
                const edge& e = meshEdges[edgeI];
                meshTools::writeOBJ(str, e.centre(meshPoints));
            }
        }

        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement(...)" << nl
            << "Writing centres of edges to split to file " << str.name()
            << endl;
    }


    // PART 3: Collect faces for refinement. Faces need to be collected in two
    // distinct categories:
    // 1. Faces found on empty patches that will be split into n faces (where n
    //    is the number of edges per face),
    // 2. Faces not on empty patches that will be always split into two
    //    faces. For each of these faces, collect the two edges found on
    //    opposing sides of the empty patch.

    if (debug)
    {
        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement" << nl
            << "Allocating face midpoints and collecting faces that are"
            << " not on empty patch."
            << endl;
    }

    // Get face anchor level based on the face type. For split face found on
    // empty patch, it is guaranteeed that we will have at least 3 points with
    // level <= anchor level. For split face not on empty patch, it is
    // guaranteed that we will have at least 4 points with level <= anchor
    // level. These are the corner points.
    labelList faceAnchorLevel(mesh_.nFaces());
    for (label faceI = 0; faceI < mesh_.nFaces(); ++faceI)
    {
        if (faceOnEmptyPatch[faceI])
        {
            // Face on empty patch, at least 3 points need to have
            // level <= anchor level
            faceAnchorLevel[faceI] = getAnchorLevel(faceI, 3);
        }
        else
        {
            // Face not on empty patch, at least 4 points need to have
            // level <= anchor level
            faceAnchorLevel[faceI] = getAnchorLevel(faceI, 4);
        }
    }

    // Split faces on empty patches will be collected in faceMidPoint list:
    // Not refined = -1
    // Shall be refined > -1 (label of added mid point)
    labelList faceMidPoint(mesh_.nFaces(), -1);

    // Split faces not on empty patches will be collected into splitFacesIntoTwo
    // dynamic list. For each of these faces, we also need to collect its two
    // edges that are found on empty patch

    // Allocate enough storage to prevent excessive resizing
    const label nSplitFacesIntoTwo = 3*cellsToRefine_.size();
    DynamicList<label> splitFacesIntoTwo(nSplitFacesIntoTwo);
    DynamicList<Pair<label> > splitFacesEmptyEdges(nSplitFacesIntoTwo);

    // Get necessary mesh data
    const labelList& meshFaceOwner = mesh_.faceOwner();
    const labelList& meshFaceNeighbour = mesh_.faceNeighbour();

    const label nFaces = mesh_.nFaces();
    const label nInternalFaces = mesh_.nInternalFaces();

    // Internal faces: look at cells on both sides. Uniquely determined since
    // the face itself is guaranteed to be same level as most refined neighbour
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        // Get owner data
        const label& own = meshFaceOwner[faceI];
        const label& ownLevel = cellLevel_[own];
        const label newOwnLevel = ownLevel + (refineCellsMask[own] ? 1 : 0);

        // Get neighbour data
        const label& nei = meshFaceNeighbour[faceI];
        const label& neiLevel = cellLevel_[nei];
        const label newNeiLevel = neiLevel + (refineCellsMask[nei] ? 1 : 0);

        if
        (
            newOwnLevel > faceAnchorLevel[faceI]
         || newNeiLevel > faceAnchorLevel[faceI]
        )
        {
            // Note: this is internal face so we don't need to check whether the
            // face is on empty patch. It can't be by definition

            // Does two things:
            // 1. Appends the face to splitFacesIntoTwo list
            // 2. Append the two edges on empty patch to splitFaceEmptyEdges
            //    list
            appendFaceSplitInfo
            (
                faceI,
                edgeOnEmptyPatch,
                edgeMidPoint,
                splitFacesIntoTwo,
                splitFacesEmptyEdges
            );
        } // End whether the face needs to be considered (split)
    } // End loop over all internal faces

    // Coupled patches handled like internal faces except now all information
    // from neighbour comes from across processor.
    // Boundary faces are more complicated since the boundary face can
    // be more refined than its owner (or neighbour for coupled patches)
    // (does not happen if refining/unrefining only, but does e.g. when
    // refinining and subsetting)

    // Memory management
    {
        // Create list for swapping boundary data
        labelList newNeiLevel(nFaces - nInternalFaces);

        forAll(newNeiLevel, i)
        {
            const label& own = meshFaceOwner[i + nInternalFaces];
            const label& ownLevel = cellLevel_[own];
            const label newOwnLevel = ownLevel + (refineCellsMask[own] ? 1 : 0);

            newNeiLevel[i] = newOwnLevel;
        }

        // Swap the list which now contains data from the other side
        syncTools::swapBoundaryFaceList(mesh_, newNeiLevel, false);

        forAll(newNeiLevel, i)
        {
            // Get face index
            const label faceI = i + nInternalFaces;

            // Get owner data (neighbour is available from before)
            const label& own = meshFaceOwner[faceI];
            const label& ownLevel = cellLevel_[own];
            const label newOwnLevel = ownLevel + (refineCellsMask[own] ? 1 : 0);

            if
            (
                newOwnLevel > faceAnchorLevel[faceI]
             || newNeiLevel[i] > faceAnchorLevel[faceI]
            )
            {
                if (faceOnEmptyPatch[faceI])
                {
                    // This face is on the empty patch and will be split into n
                    // faces (n is the number of edges for this face) and the
                    // face mid point will be added. Mark the face for splitting
                    faceMidPoint[faceI] = 12345;
                }
                else
                {
                    // Does two things:
                    // 1. Appends the face to splitFacesIntoTwo list
                    // 2. Append the two edges on empty patch to splitFaceEmptyEdges
                    //    list
                    appendFaceSplitInfo
                    (
                        faceI,
                        edgeOnEmptyPatch,
                        edgeMidPoint,
                        splitFacesIntoTwo,
                        splitFacesEmptyEdges
                    );
                } // End if the face is not on empty patch
            } // End whether the face needs to be considered
        } // End loop over all boundary faces
    } // End memory management for syncing owner/neighbour face levels


    // Add face points. Note: no need to sync face mid points (as we did for
    // edge mid points) since processor faces do not introduce new points, only
    // faces on empty patch do

    // Get face centres
    const vectorField& meshFaceCentres = mesh_.faceCentres();

    // Loop through faces on empty patches only
    forAll (boundaryMesh, patchI)
    {
        // Get current patch
        const polyPatch& curPatch = boundaryMesh[patchI];

        // Check whether this patch is emptyPolyPatch
        if (isA<emptyPolyPatch>(curPatch))
        {
            // Get start and face labels
            const label startFaceI = curPatch.start();
            const label endFaceI = startFaceI + curPatch.size();

            // Loop through all empty patch faces (global indexing)
            for (label faceI = startFaceI; faceI < endFaceI; ++faceI)
            {
                if (faceMidPoint[faceI] > -1)
                {
                    // Face on empty patch marked to be split. Add the point at
                    // face centre and replace faceMidPoint with new point label

                    faceMidPoint[faceI] = ref.setAction
                    (
                        polyAddPoint
                        (
                            meshFaceCentres[faceI], // Point
                            -1,                     // No master ID
                            -1,                     // Zone for point
                            true                    // Supports a cell
                        )
                    );
                } // End if face marked for splitting
            } // End loop over all faces on empty patch
        } // End if emptyPolyPatch check
    } // End loop for all patches

    // Write out all split faces as a face set for debugging
    if (debug)
    {
        // Create faceSet containing all faces that need to be split into n
        // faces (n is the number of edges on the face)
        faceSet splitNFaces(mesh_, "splitNFaces", 3*cellsToRefine_.size());

        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] > -1)
            {
                splitNFaces.insert(faceI);
            }
        }

        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement(...)" << nl
            << "Writing " << splitNFaces.size()
            << " faces to split in N to faceSet " << splitNFaces.objectPath()
            << endl;

        splitNFaces.write();

        // Create faceSet containing all faces that need to be split into 2
        // faces (faces not on empty patch)
        faceSet splitTwoFaces(mesh_, "splitTwoFaces", 3*cellsToRefine_.size());

        forAll (splitFacesIntoTwo, i)
        {
            // Insert face index into splitTwoFaces
            splitTwoFaces.insert(splitFacesIntoTwo[i]);
        }

        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement(...)" << nl
            << "Writing " << splitTwoFaces.size()
            << " faces to split in two to faceSet "
            << splitTwoFaces.objectPath() << endl;

        splitTwoFaces.write();
    }


    // Now we have all the information we need to perform the refinement and we
    // no longer need to refer to cellsToRefine_. The information is in:
    // - refineCellsMask = true : cell needs to be split
    // - edgeMidPoint >= 0     : edge on empty patch that needs to be split
    // - faceMidPoint >= 0     : face on empty patch that needs to be split
    //                           into n faces (where n is the number of edges)
    // - splitFacesIntoTwo     : list of faces that need to be split into two
    //                           (face not on empty patch)
    // - splitFacesEmptyEdges  : holds the two edges of the face which needs to
    //                           be split into two


    // PART 4: Get corner and anchor points for all cells

    if (debug)
    {
        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement(...)" << nl
            << "Finding cell anchorPoints" << endl;
    }

    // Get anchor points for each cell: points that have the same or lower
    // refinement level as the cell
    List<dynamicLabelList> cellAnchorPointsDynamic(mesh_.nCells());

    // Loop through all cells
    forAll(refineCellsMask, cellI)
    {
        if (refineCellsMask[cellI])
        {
            // The cell will be refined, set capacity to 8 to prevent excessive
            // resizing
            cellAnchorPointsDynamic[cellI].setCapacity(8);
        }
    }

    // Get necessary mesh data
    const labelListList& meshPointCells = mesh_.pointCells();

    // Loop through all points
    forAll(pointLevel_, pointI)
    {
        // Get point cells
        const labelList& pCells = meshPointCells[pointI];

        // Loop through all cells sharing this point
        forAll(pCells, pCellI)
        {
            // Get current cell index
            const label& cellI = pCells[pCellI];

            if
            (
                refineCellsMask[cellI]
             && pointLevel_[pointI] <= cellLevel_[cellI]
            )
            {
                // This point cell is marked for refinement and its point level
                // is smaller or equal to cell level, append the point
                cellAnchorPointsDynamic[cellI].append(pointI);
            }
        }
    }

    // Loop through all cells and check whether at least 6 anchor points
    // have been found (minimum requirement for a triangular prism)

    // Collect cellAnchorPoint into a List<labelList> instead of
    // List<dynamicList>
    labelListList cellAnchorPoints(mesh_.nCells());

    // Get cell points for error output
    const labelListList& meshCellPoints = mesh_.cellPoints();

    forAll(refineCellsMask, cellI)
    {
        // First some sanity checks
        if (refineCellsMask[cellI])
        {
            // Cell selected for refinement
            if (cellAnchorPointsDynamic[cellI].size() < 6)
            {
                // Cell has less than 6 anchor points. Issue an error and report
                // cell points
                const labelList& cPoints = meshCellPoints[cellI];

                FatalErrorIn
                (
                    "prismatic2DRefinement::setPrismatic2DRefinement(...)"
                )   << "Cell " << cellI
                    << " of level " << cellLevel_[cellI]
                    << " does not seem to have enough points of "
                    << " lower level" << endl
                    << "cellPoints:" << cPoints << endl
                    << "pointLevels:"
                    << IndirectList<label>(pointLevel_, cPoints)() << endl
                    << abort(FatalError);
            }
            else if (cellAnchorPointsDynamic[cellI].size() % 2 != 0)
            {
                // Cell has odd number of anchor points. This is not allowed and
                // indicates an invalid mesh
                const labelList& cPoints = meshCellPoints[cellI];

                FatalErrorIn
                (
                    "prismatic2DRefinement::setPrismatic2DRefinement(...)"
                )   << "Cell " << cellI
                    << " of level " << cellLevel_[cellI]
                    << " has odd number of anchor points"
                    << " (should be even for 2D mesh). "
                    << "cellPoints:" << cPoints << endl
                    << "pointLevels:"
                    << IndirectList<label>(pointLevel_, cPoints)() << endl
                    << abort(FatalError);
            }
        }

        // Tranfer the dynamic list for each cell into an ordinary list
        cellAnchorPoints[cellI].transfer(cellAnchorPointsDynamic[cellI]);
    }


    // PART 5: Add the cells

    if (debug)
    {
        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement(...)" << nl
            << " Adding cells."
            << endl;
    }

    // We should have exactly n new cells per each split cell, where n is the
    // number of anchor points in a cell divided by two. In order to determine
    // owner/neighbours of new and modified faces, we need to know which cell
    // came from which point. The mapping is not uniquely defined as in
    // polyhedralRefinement when we had 1 point = 1 cell. Here, we have two
    // points that correspond to a single cell, one on one side of the empty
    // patch and the other on other side. This information will be collected in
    // a HashTable<label, Pair<label> >, where the key will be a pair of
    // global point index and global cell index, while the value is local index
    // into cellAddedCells list
    labelListList cellAddedCells(mesh_.nCells());
    HashTable<label, Pair<label>, Hash<FixedList<label, 2> > >
        pointCellToAddedCellMap(6*cellsToRefine_.size());

    // Get mesh data
    const cellZoneMesh& cellZones = mesh_.cellZones();
    const cellList& meshCells = mesh_.cells();
    const faceList& meshFaces = mesh_.faces();
    const labelListList& meshPointEdges = mesh_.pointEdges();


    // Loop through all faces
    forAll(faceMidPoint, faceI)
    {
        // Check whether this face needs to be split into n faces
        if (faceMidPoint[faceI] > -1)
        {
            // This is a face that will be split and is on empty patch by
            // definition, get the cell index by looking at owner only
            const label& cellI = meshFaceOwner[faceI];

            // Get cell added cells
            labelList& cAdded = cellAddedCells[cellI];

            // Check whether the cell added cells are empty. This means that we
            // haven't visited the first face yet. If it is not empty, we have
            // already visited one face, which is enough
            if (cAdded.empty())
            {
                // First face that hasn't been visited. Start adding cells
                // point-by-poing and keep track of mapping necessary for
                // splitting other (not on empty patch) faces into two

                // Set the total number of added cells to number of anchors
                // divided by two. Note: number of anchors needs to be an even
                // number (6 for triangular prism, 8 for hex, 10 for pentagonal
                // prism, etc.)
                const labelList& cAnchors = cellAnchorPoints[cellI];
                cAdded.setSize(cAnchors.size()/2);

                // Helper variable to distinguish between first and successive
                // cells (first will have the original index)
                label cellCounter = 0;

                // Get current face
                const face& f = meshFaces[faceI];

                // Loop through face points
                forAll (f, fpI)
                {
                    // Get point index
                    const label& pointI = f[fpI];

                    // Find anchor point in local list if present
                    const label anchorI = findIndex(cAnchors, pointI);

                    if (anchorI != -1)
                    {
                        // This point is anchor, add the cell

                        if (cellCounter == 0)
                        {
                            // This is first cell, simply set the existing index
                            cAdded[cellCounter] = cellI;

                            // Update refinement level indicator field to 1
                            // since the original cell will be refined
                            refinementLevelIndicator_[cellI] = REFINED;
                        }
                        else
                        {
                            // Other cells, need to add the cells
                            cAdded[cellCounter] = ref.setAction
                            (
                                polyAddCell
                                (
                                    -1,                         // M. point
                                    -1,                         // M. edge
                                    -1,                         // M. face
                                    cellI,                      // M. cell
                                    cellZones.whichZone(cellI)  // M. zone
                                )
                            );
                        }

                        // Collect the point-cell mapping into local index
                        // of cell added cells for point on this side
                        pointCellToAddedCellMap.insert
                        (
                            Pair<label>(pointI, cellI),
                            cellCounter
                        );

                        // This is only one side, we need to also collect
                        // the other point on the other side. Get point edges
                        // for this point
                        const labelList& pEdges =
                            meshPointEdges[pointI];

                        // Loop through point edges
                        forAll (pEdges, peI)
                        {
                            // Get the edge index
                            const label& edgeI = pEdges[peI];

                            if (!edgeOnEmptyPatch[edgeI])
                            {
                                // Edge is not on empty patch, therefore
                                // this is the edge we're looking
                                // for. Collect the other point of the edge
                                const label pointJ =
                                    meshEdges[edgeI].otherVertex(pointI);

                                // Collect the point-cell mapping into local
                                // index of cell added cells for this point
                                pointCellToAddedCellMap.insert
                                (
                                    Pair<label>(pointJ, cellI),
                                    cellCounter
                                );
                            }
                        }

                        // Now we have finished adding the cell and also
                        // adding the necessary mapping for this added cell

                        // Increment cell counter
                        ++cellCounter;

                    } // Else point is not anchor: nothing to do
                } // End loop over all face points

                // Sanity check: number of counted cells must be
                // equal to size of cellAddedCells. This means that
                // we have correctly marked the anchor points
                if (cellCounter != cAdded.size())
                {
                    FatalErrorIn
                    (
                        "prismatic2DRefinement::setPrismatic2DRefinement(...)"
                    )   << "Problem while adding cells."
                        << nl
                        << "Going through base face on empty patch and adding"
                        << " cells, we collected: " << cellCounter << " cells."
                        << nl
                        << "While the number of anchor points is: "
                        << cAnchors.size()
                        << nl
                        << "The number of added cells based on number of anchor"
                        << " points is: "
                        << cAdded.size()
                        << nl
                        << "Additional information: "
                        << nl
                        << "cellI: " << cellI << ", faceI: " << faceI
                        << abort(FatalError);
                }

            } // End if cell added cells empty
        } // End if face needs to be split into n
    } // End loop over all faces


    // PART 6: Adding faces

    // 6.1. Existing faces on empty patches that get split (into n faces
    //      where n is the number of points or edges)
    // 6.2. Existing faces not on empty patches that get split into two
    // 6.3. Existing faces that do not get split but only edges get split
    // 6.4. Existing faces that do not get split but get new owner/neighbour
    // 6.5. New internal faces inside split cells

    if (debug)
    {
        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement(...)" << nl
            << " Marking faces to be handled"
            << endl;
    }

    // Get all faces to split:
    // a) All faces of a cell being split
    // b) All faces on empty patch that are being split
    // c) All faces not on empty patch that are being split
    // d) Both faces of an edge that is being split
    // Note: although a bit redundant, loop over everything above
    boolList facesToSplit(mesh_.nFaces(), false);

    // Also collect all split faces which will be needed in 6.3
    boolList allSplitFaces(mesh_.nFaces(), false);

    // Get edge faces
    const labelListList& meshEdgeFaces = mesh_.edgeFaces();

    // a) All faces of a cell that is being split
    forAll(refineCellsMask, cellI)
    {
        if (refineCellsMask[cellI])
        {
            const cell& cFaces = meshCells[cellI];

            forAll(cFaces, i)
            {
                facesToSplit[cFaces[i]] = true;
            }
        }
    }

    // b) All faces on empty patch that are being split
    forAll(faceMidPoint, faceI)
    {
        if (faceMidPoint[faceI] > -1)
        {
            // Mark face in both lists
            facesToSplit[faceI] = true; // Used through 6.1-6.5
            allSplitFaces[faceI] = true; // Used in 6.3
        }
    }

    // c) All faces not on empty patch that are being split
    forAll(splitFacesIntoTwo, i)
    {
        // Get face index
        const label& faceI = splitFacesIntoTwo[i];

        // Mark face in both lists
        facesToSplit[faceI] = true;
        allSplitFaces[faceI] = true;
    }

    // d) Both faces of an edge that is being split
    forAll(edgeMidPoint, edgeI)
    {
        if (edgeMidPoint[edgeI] > -1)
        {
            const labelList& eFaces = meshEdgeFaces[edgeI];

            forAll(eFaces, i)
            {
                facesToSplit[eFaces[i]] = true;
            }
        }
    }

    // Note: after splitting a certain face during parts 6.1. to 6.4.,
    // facesToSplit for that face will be set back to 0, i.e. marked as finished


    // PART 6.1. Add/modify faces for each face on empty patch that is being
    // split

    if (debug)
    {
        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement(...)" << endl
            << "Splitting faces on empty patches..." << endl;
    }

    forAll(faceMidPoint, faceI)
    {
        if (faceMidPoint[faceI] > -1 && facesToSplit[faceI])
        {
            // Face has not been split.
            // Note: although facesToSplit can't be different than 1 here and
            // the second check can be ommitted, it is left for clarity

            // Get the face
            const face& f = meshFaces[faceI];

            // Flag to control whether the original faceI has been used
            // Note: original face gets modified, other n - 1 faces are added,
            // where n is the number of points/edges of a face
            bool modifiedFace = false;

            // Get anchor level for the face
            const label& anchorLevel = faceAnchorLevel[faceI];

            // New face always has four points/edges for arbitrary polygon
            face newFace(4);

            // Loop through all points of original face
            forAll(f, fpI)
            {
                const label& pointI = f[fpI];

                if (pointLevel_[pointI] <= anchorLevel)
                {
                    // This point is anchor, start collecting face

                    // Create dynamic list (because of append) for face vertices
                    // and append the first (anchor) point
                    dynamicLabelList faceVerts(4);
                    faceVerts.append(pointI);

                    // Walk forward to mid point.
                    // - if next is +2 midpoint is +1
                    // - if next is +1 it is midpoint
                    // - if next is +0 there has to be edgeMidPoint

                    // Appends all points from this point to face mid point
                    walkFaceToMid
                    (
                        edgeMidPoint,
                        anchorLevel,
                        faceI,
                        fpI,
                        faceVerts
                    );

                    // Append face mid point
                    faceVerts.append(faceMidPoint[faceI]);

                    // Append all points from face mid point to starting point
                    walkFaceFromMid
                    (
                        edgeMidPoint,
                        anchorLevel,
                        faceI,
                        fpI,
                        faceVerts
                    );

                    // Transfer dynamic list to a face (ordinary list)
                    newFace.transfer(faceVerts);
                    faceVerts.clear();

                    // Set new owner/neighbour indices based on split cells
                    label own, nei;
                    setNewFaceNeighbours
                    (
                        pointCellToAddedCellMap,
                        cellAddedCells,
                        faceI,
                        pointI, // Anchor point index

                        own,
                        nei
                    );

                    if (debug)
                    {
                        // Check orientation of the split face for debugging
                        checkNewFaceOrientation(ref, faceI, newFace);
                    }


                    // Finally insert the modification/addition instruction into
                    // the topo changer engine
                    if (!modifiedFace)
                    {
                        // Modify first face
                        modifiedFace = true;
                        modifyFace(ref, faceI, newFace, own, nei);
                    }
                    else
                    {
                        // Add additional faces
                        addFace(ref, faceI, newFace, own, nei);
                    }
                } // End point anchor check
            } // End for all points

            // Mark face as handled
            facesToSplit[faceI] = false;

        } // End if this face needs to be split
    } // End for all faces


    // PART 6.2. Add/modify faces for each face not on empty patch that is being
    // split into two

    if (debug)
    {
        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement(...)" << nl
            << "Splitting faces not on empty patches..."
            << endl;
    }

    // Loop through faces that are not on empty patch. These will be split into
    // two faces only
    forAll (splitFacesIntoTwo, i)
    {
        // Get face index
        const label& faceI = splitFacesIntoTwo[i];

        // Check whether the face is marked for splitting. A bit redundant but
        // will be left for clarity
        if (facesToSplit[faceI])
        {
            // Face has not been split, grab the face
            const face& f = meshFaces[faceI];

            // Additional sanity check
            if (f.size() != 4)
            {
                FatalErrorIn
                (
                    "prismatic2DRefinement::appendFaceSplitInfo(...)"
                )   << "The original face has: " << f.size() << " points,"
                    << " while it should have exactly 4 points in order"
                    << " to split it in two."
                    << " faceI: " << faceI
                    << abort(FatalError);
            }

            // Flag to control whether the original faceI has been used
            // Note: original face gets modified, other face gets added
            bool modifiedFace = false;

            // Get anchor level for the face
            const label& anchorLevel = faceAnchorLevel[faceI];

            // Mark visited points to avoid adding the face twice
            boolList visitedPoint(f.size(), false);

            // New face always has four points/edges for 2D face splitting of
            // a face that is not on empty patch
            face newFace(4);

            // Loop through all points of original face
            forAll (f, fpI)
            {
                // Get point index
                const label& pointI = f[fpI];

                if
                (
                    !visitedPoint[fpI]
                 && pointLevel_[pointI] <= anchorLevel
                )
                {
                    // This point is anchor and it hasn't been visited yet,
                    // start collecting face

                    // Collect the new face in the following order:
                    // 1. This point
                    // 2. Edge mid points for the edge that contains this point
                    //    and the edge that is on the other side
                    // 3. Other point (on the other side)

                    // 1. Set this point and mark it as visited
                    newFace[0] = pointI;
                    visitedPoint[fpI] = true;

                    // 2. Get edge mid point for edge containing this point
                    // Fetch the two edges on both sides
                    const Pair<label>& edgesOnOppositeSides =
                        splitFacesEmptyEdges[i];

                    // Get the edge indices and edges
                    const label& edgeIndexI = edgesOnOppositeSides.first();
                    const label& edgeIndexJ = edgesOnOppositeSides.second();

                    const edge& edgeI = meshEdges[edgeIndexI];
                    const edge& edgeJ = meshEdges[edgeIndexJ];

                    // Additional sanity check
                    if
                    (
                        (edgeMidPoint[edgeIndexI] == -1)
                     || (edgeMidPoint[edgeIndexJ] == -1)
                    )
                    {
                        // Edges are not marked for refinement, issue an error
                        FatalErrorIn
                        (
                            "prismatic2DRefinement::appendFaceSplitInfo(...)"
                        )   << "Trying to split a face into two, but"
                            << " edges on empty patches are not properly set."
                            << nl
                            << "Edge: " << edgeIndexI << " with new point: "
                            << edgeMidPoint[edgeIndexI]
                            << nl
                            << "Edge: " << edgeIndexJ << " with new point: "
                            << edgeMidPoint[edgeIndexJ]
                            << abort(FatalError);
                    }

                    if ((edgeI.start() == pointI) || (edgeI.end() == pointI))
                    {
                        // Current point is on edgeI, set edgeI midpoint and
                        // then edgeJ midpoint
                        newFace[1] = edgeMidPoint[edgeIndexI];
                        newFace[2] = edgeMidPoint[edgeIndexJ];
                    }
                    else if
                    (
                        (edgeJ.start() == pointI) || (edgeJ.end() == pointI)
                    )
                    {
                        // Current is on edgeJ, set edgeJ midpoint and then
                        // edgeI midpoint
                        newFace[1] = edgeMidPoint[edgeIndexJ];
                        newFace[2] = edgeMidPoint[edgeIndexI];
                    }
                    else
                    {
                        // Point not on either of edges, issue an error
                        FatalErrorIn
                        (
                            "prismatic2DRefinement::appendFaceSplitInfo(...)"
                        )   << "Trying to split a face into two, but"
                            << " the point: " << pointI << " can't be found"
                            << " on either of edges. "
                            << nl
                            << "Edge: " << edgeIndexI << " with new point: "
                            << edgeMidPoint[edgeIndexI]
                            << nl
                            << "Edge: " << edgeIndexJ << " with new point: "
                            << edgeMidPoint[edgeIndexJ]
                            << abort(FatalError);
                    }

                    // At this point, we have added three points: original
                    // point, first edge mid point and second edge mid point.

                    // 3. Add the other point
                    // Get point edges for this point
                    const labelList& pEdges = meshPointEdges[pointI];

                    // Loop through all edges
                    forAll (pEdges, peI)
                    {
                        // Get the edge index
                        const label& pointEdgeI = pEdges[peI];

                        if (!edgeOnEmptyPatch[pointEdgeI])
                        {
                            // Edge is not on empty patch, therefore this is the
                            // edge we're looking for. Collect the other point
                            const label pointJ =
                                meshEdges[pointEdgeI].otherVertex(pointI);

                            // Insert the point into the face at the last location
                            newFace[3] = pointJ;

                            // Mask local point index as visited by going
                            // through the face again
                            forAll (f, fpJ)
                            {
                                if (f[fpJ] == pointJ)
                                {
                                    // Found local index of the point, mask it
                                    visitedPoint[fpJ] = true;
                                }
                            }
                        }
                    }

                    // The face is now complete, set new owner/neighbour indices
                    // based on split cells
                    label own, nei;

                    // Set new face owner/neighbour pair
                    setNewFaceNeighbours
                    (
                        pointCellToAddedCellMap,
                        cellAddedCells,
                        faceI,
                        pointI, // Anchor point index

                        own,
                        nei
                    );

                    // We need to revert the face if the edge between this point
                    // and the next point is not split. This follows from
                    // definition of face as ordered set of points (defining
                    // orientation) and the splitting procedure. Note: edge
                    // ordering in face is the same as point ordering so the
                    // point index can be used as first face edge index
                    if (edgeMidPoint[meshFaceEdges[faceI][fpI]] == -1)
                    {
                        newFace = newFace.reverseFace();
                    }

                    if (debug)
                    {
                        // Check orientation of the split face for debugging
                        checkNewFaceOrientation(ref, faceI, newFace);
                    }


                    // Finally insert the modification/addition instruction into
                    // the topo changer engine
                    if (!modifiedFace)
                    {
                        // Modify first face
                        modifiedFace = true;
                        modifyFace(ref, faceI, newFace, own, nei);
                    }
                    else
                    {
                        // Add additional faces
                        addFace(ref, faceI, newFace, own, nei);
                    }
                } // End if point is anchored and has not been visited
            } // End loop over all face points

            // Mark face as handled
            facesToSplit[faceI] = false;

        } // End if face is split
    } // End for all faces that should be split into two


    // PART 6.3. Modify faces that do not get split but have edges that are
    // being split

    if (debug)
    {
        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement(...)" << nl
            << "Modifying faces with split edges..."
            << endl;
    }

    forAll(edgeMidPoint, edgeI)
    {
        if (edgeMidPoint[edgeI] > -1)
        {
            // This is an edge that is going to be split, get edge faces
            const labelList& eFaces = meshEdgeFaces[edgeI];

            // Loop through all faces of an edge
            forAll(eFaces, i)
            {
                // Get face index
                const label faceI = eFaces[i];

                // Check whether this is not a face that's been split and that
                // the face has not been handled yet. The second check is
                // necessary since we go through edge faces instead of just
                // faces
                if (!allSplitFaces[faceI] && facesToSplit[faceI])
                {
                    // This is unsplit face that has not been handled

                    // Get face and face edges
                    const face& f = meshFaces[faceI];
                    const labelList& fEdges = meshFaceEdges[faceI];

                    // Create a dynamic list containing new face vertices
                    dynamicLabelList newFaceVerts(f.size());

                    // Append all original points and all edge mid points
                    forAll(f, fpI)
                    {
                        newFaceVerts.append(f[fpI]);

                        const label edgeI = fEdges[fpI];

                        if (edgeMidPoint[edgeI] > -1)
                        {
                            newFaceVerts.append(edgeMidPoint[edgeI]);
                        }
                    }

                    // Create a face from dynamic list by transfer
                    face newFace(newFaceVerts.xfer());


                    // The point with the lowest level should be an anchor
                    // point of the neighbouring cells.
                    const label anchorFpI = findMinLevel(f);

                    label own, nei;
                    setNewFaceNeighbours
                    (
                        pointCellToAddedCellMap,
                        cellAddedCells,
                        faceI,
                        f[anchorFpI], // Anchor point index

                        own,
                        nei
                    );


                    if (debug)
                    {
                        // Check orientation of the new face for debugging
                        checkNewFaceOrientation(ref, faceI, newFace);
                    }

                    // Modify the face
                    modifyFace(ref, faceI, newFace, own, nei);

                    // Mark face as handled
                    facesToSplit[faceI] = false;

                } // End if unsplit, unhandled face
            } // End for all edge faces
        } // End if edge has been cut
    } // End for all edges


    // PART 6.4: Modify faces that do not get split but whose owner/neighbour
    // change due to splitting

    if (debug)
    {
        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement(...)" << nl
            << " Changing owner/neighbour for otherwise unaffected faces..."
            << endl;
    }

    forAll(facesToSplit, faceI)
    {
        // All remaining unnaffected faces are the ones whose owner/neighbour
        // changed
        if (facesToSplit[faceI])
        {
            // Get the face
            const face& f = meshFaces[faceI];

            // The point with the lowest level should be an anchor point of the
            // neighbouring cells
            label anchorFpI = findMinLevel(f);

            label own, nei;
            setNewFaceNeighbours
            (
                pointCellToAddedCellMap,
                cellAddedCells,
                faceI,
                f[anchorFpI], // Anchor point

                own,
                nei
            );

            // Modify the face, changing owner and neighbour
            modifyFace(ref, faceI, f, own, nei);

            // Mark face as handled
            facesToSplit[faceI] = false;

        } // End if the face needs to be handled
    } // End for all faces


    // PART 6.5. Add new internal faces inside split cells

    if (debug)
    {
        Pout<< "prismatic2DRefinement::setPrismatic2DRefinement(...)" << nl
            << " Adding new internal faces for split cells..."
            << endl;
    }

    // Mark-up filed for visited cells (since we are going through faces)
    boolList cellsToSplit(mesh_.nCells(), true);

    // Loop through faces in the same way as we did when we were adding
    // cells. This order is important since it ensures easy determination of
    // owner/neighbour cells for new faces
    forAll(faceMidPoint, faceI)
    {
        // Get owner of the face. For face on empty patch
        const label& cellI = meshFaceOwner[faceI];

        // Check whether this face has been split and whether the cell has been
        // handled (internal faces already created for this cell)
        if (faceMidPoint[faceI] > -1 && cellsToSplit[cellI])
        {
            // Face is split and hasn't been visited yet. Get the face and edges
            const face& f = meshFaces[faceI];
            const labelList& fEdges = meshFaceEdges[faceI];

            // Get anchor points for this cell and cell added cells
            const labelList& cAnchors = cellAnchorPoints[cellI];
            const labelList& cAdded = cellAddedCells[cellI];

            // Count number of added faces (helper variable to determine
            // owner/neighbour)
            label nAddedFaces = 0;

            // Loop through face points
            forAll (f, fpI)
            {
                // Get point index
                const label& pointI = f[fpI];

                // If this point is not an anchor, it has already been handled
                // (by going through anchors), skip
                if (findIndex(cAnchors, pointI) == -1)
                {
                    continue;
                }

                // Get corresponding edge index and edge (between fpI and
                // fpI + 1) by definition of faceEdges
                const label& edgeI = fEdges[fpI];
                const edge& e = meshEdges[edgeI];

                // Grab other point
                const label& pointJ = e.otherVertex(pointI);

                if (pointJ == -1)
                {
                    // If pointJ is equal to -1, this means that the pointI
                    // was not found on edge, something went wrong
                    FatalErrorIn
                    (
                        "prismatic2DRefinement::setPrismatic2DRefinement(...)"
                    )   << "Point: " << pointI << " not found on edge: "
                        << edgeI << nl
                        << "Looping through face points and face edges did"
                        << " not ensure synchronous behaviour."
                        << abort(FatalError);
                }

                // Create the new face
                face newFace(4);

                // Note: there are three possible variants:
                //   i) Edge is split and the other point is anchor. Collection
                //      of the face starts from edgeMidPoint
                //  ii) Edge is split and the other point is not an
                //      anchor. Collection of the face starts from other point
                // iii) Edge is not split and the other point is not an
                //      anchor. Collection of the face starts from other point
                // Variants ii) and iii) can be handled together, while variant
                // i) has to be handled separately.

                // Whether the edge is split
                const bool isEdgeSplit = edgeMidPoint[edgeI] > -1;

                // Whether the other point is anchor or not
                const bool isOtherEdgePointAnchor
                    = findIndex(cAnchors, pointJ) > -1;

                // Check if the edge is split and whether the other edge point
                // is an anchor
                if (isEdgeSplit && isOtherEdgePointAnchor)
                {
                    // Variant i) Edge is split and other edge point is an anchor

                    // Create the new face and start collecting points
                    // a) edgeMidPoint for this edge
                    // b) faceMidPoint for this face
                    // c) faceMidPoint for the face on the other side
                    // d) edgeMidPoint for the edge on the other side

                    // a) edgeMidPoint for this edge
                    newFace[0] = edgeMidPoint[edgeI];

                    // b) and c): adding both face mids
                    addFaceMids
                    (
                        faceMidPoint,
                        faceOnEmptyPatch,
                        faceI,
                        cellI,
                        newFace
                    );

                    // d) edgeMidPoint for the edge on the other side
                    // The other edge is uniquely defined as the edge on empty
                    // patch sharing the same face as this edge

                    // Get the edge faces
                    const labelList& eFaces = meshEdgeFaces[edgeI];

                    // Loop through edge faces
                    forAll(eFaces, i)
                    {
                        // Get the face and check whether it is on empty patch
                        const label& faceK = eFaces[i];

                        if (!faceOnEmptyPatch[faceK])
                        {
                            // Found the face, need to search its edges
                            const labelList& otherFaceEdges =
                                meshFaceEdges[faceK];

                            forAll(otherFaceEdges, j)
                            {
                                // Get the edge
                                const label& edgeJ = otherFaceEdges[j];

                                if (edgeOnEmptyPatch[edgeJ] && (edgeI != edgeJ))
                                {
                                    // Edge is on empty patch, this must be the
                                    // one we are looking for. Add its midpoint
                                    // and double check if it is valid
                                    if (edgeMidPoint[edgeJ] > -1)
                                    {
                                        newFace[3] = edgeMidPoint[edgeJ];
                                        break;
                                    }
                                    else
                                    {
                                        FatalErrorIn
                                        (
                                            "prismatic2DRefinement::"
                                            "setPrismatic2DRefinement(...)"
                                        )   << "Other edge: "
                                            << edgeJ
                                            << " has not been selected for splitting,"
                                            << " while the edge on original side: "
                                            << edgeI
                                            << " has been selected."
                                            << abort(FatalError);
                                    }
                                } // End if this is our "other" edge
                            } // End for all other (non empty patch) face edges

                            // Break out since we must have found the candidate
                            break;

                        } // End if face not on empty patch

                    } // End for all edge faces

                } // End if this edge is split and the other point is anchor
                else if (!isOtherEdgePointAnchor)
                {
                    // Variants ii) and iii). Either the edge is split and the
                    // other point is not an anchor or the edge is not split and
                    // the other point is not an anchor

                    // Create the new face and start collecting points
                    // a) other point of this edge
                    // b) faceMidPoint for this face
                    // c) faceMidPoint for the face on the other side
                    // d) other point on the other side

                    // a) other point of this edge
                    newFace[0] = pointJ;

                    // b) and c): adding both face mids
                    addFaceMids
                    (
                        faceMidPoint,
                        faceOnEmptyPatch,
                        faceI,
                        cellI,
                        newFace
                    );

                    // d) other point on the other side
                    // The other point is uniquely defined as the other point of
                    // the edge of this point which is not on empty patch

                    // Get point edges
                    const labelList& pEdges = meshPointEdges[pointJ];

                    // Loop through all edges
                    forAll(pEdges, i)
                    {
                        // Get the edge index
                        const label& edgeJ = pEdges[i];

                        if (!edgeOnEmptyPatch[edgeJ])
                        {
                            // Found our edge, set the point on the other side
                            // of the edge as the last point in face
                            newFace[3] = meshEdges[edgeJ].otherVertex(pointJ);
                            break;
                        }
                    } // End loop over all point edges

                } // End if the other point is not an anchor
                else
                {
                    // The edge is not split and the other point is an
                    // anchor. This should never happen
                    FatalErrorIn
                    (
                        "prismatic2DRefinement::setPrismatic2DRefinement(...)"
                    )   << "Attempted to create internal face for an edge that"
                        << " is not split and the other point that is an anchor."
                        << nl
                        << "Cell: " << cellI
                        << ", point: " << pointI
                        << ", other edge point: " << pointJ
                        << nl
                        << "Anchor points for cell are: " << cAnchors
                        << abort(FatalError);
                }

                // Now we have the face defined, set owner and neighbour.
                // Note: owner and neighbour are uniquely defined since we have
                // gone through the face in the same way as we did while adding
                // cells. This ensured easy definition of owner/neighbour cells
                const label own = cAdded[nAddedFaces];
                const label nei =
                    nAddedFaces < cAdded.size() - 1
                  ? cAdded[nAddedFaces + 1]
                  : cAdded[0];

                // According to the definition of adding faces, the first n - 1
                // faces need to be reverted, while the last one is correctly
                // oriented
                if (nAddedFaces < cAdded.size() - 1)
                {
                    newFace = newFace.reverseFace();
                }


                // Debug: check orientation
                if (debug)
                {
                    // Get owner/neighbour points
                    point ownPt, neiPt;

                    if (nAddedFaces < cAdded.size() - 1)
                    {
                        // Original owner/neighbour
                        ownPt = meshPoints[pointI];
                        neiPt = meshPoints[pointJ];
                    }
                    else
                    {
                        // Flipped owner/neighbour for last face
                        ownPt = meshPoints[pointJ];
                        neiPt = meshPoints[pointI];
                    }

                    checkInternalOrientation
                    (
                        ref,
                        cellI,
                        faceI,
                        ownPt,
                        neiPt,
                        newFace
                    );
                }

                // Finally, add the face. Note: ignoring return of new face
                // index from ref.setAction(polyAddFace(...)) call
                ref.setAction
                (
                    polyAddFace
                    (
                        newFace, // face
                        own,     // owner
                        nei,     // neighbour
                        -1,      // master point
                        -1,      // master edge
                        0,       // master face for addition
                        false,   // flux flip
                        -1,      // patch for face
                        -1,      // zone for face
                        false    // face zone flip
                    )
                );

                // Increment number of added faces
                ++nAddedFaces;

            } // End loop over all point (and edges) of the face

            // Finished adding internal faces. Mark the cell as handled
            cellsToSplit[cellI] = false;

        } // End if face is split into n and cell has not been handled
    } // End for all faces

    // Debug: check minimum point index of added points, needs to be equal to
    // number of points in the original mesh
    if (debug)
    {
        label minPointI = labelMax;
        label maxPointI = labelMin;

        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] > -1)
            {
                minPointI = min(minPointI, faceMidPoint[faceI]);
                maxPointI = max(maxPointI, faceMidPoint[faceI]);
            }
        }
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] > -1)
            {
                minPointI = min(minPointI, edgeMidPoint[edgeI]);
                maxPointI = max(maxPointI, edgeMidPoint[edgeI]);
            }
        }

        if (minPointI != labelMax && minPointI != mesh_.nPoints())
        {
            FatalErrorIn("prismatic2DRefinement::setPrismatic2DRefinement(...)")
                << "Added point labels not consecutive to existing mesh points."
                << nl
                << "mesh_.nPoints():" << mesh_.nPoints()
                << " minPointI: " << minPointI
                << " maxPointI: " << maxPointI
                << abort(FatalError);
        }
    }
}


void Foam::prismatic2DRefinement::setPrismatic2DUnrefinement
(
    polyTopoChange& ref
) const
{
    // Note: assumes that splitPointsToUnrefine_ are set prior to the function
    // call

    // Check whether the refinementLevelIndicator is valid
    if (refinementLevelIndicator_.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "prismatic2DRefinement::setPrismatic2DUnrefinement(...)"
        )   << "Refinement level indicator list has invalid size: "
            << refinementLevelIndicator_.size()
            << ", number of cells: " << mesh_.nCells()
            << nl
            << "Make sure to call setPrismatic2DRefinement(...) before"
            << " calling setPrismatic2DUnrefinement(...)."
            << abort(FatalError);
    }

    // Get point cells necessary for debug and face removal
    const labelListList& meshPointCells = mesh_.pointCells();

    if (debug)
    {
        Pout<< "prismatic2DRefinement::setPrismatic2DUnrefinement"
            << "(polyTopoChange& ref)"
            << nl
            << "Checking validity of cellLevel before setting unrefinement."
            << endl;

        forAll(cellLevel_, cellI)
        {
            if (cellLevel_[cellI] < 0)
            {
                FatalErrorIn
                (
                    "prismatic2DRefinement::setPrismatic2DUnrefinement"
                    "(polyTopoChange& ref)"
                )   << "Illegal cell level " << cellLevel_[cellI]
                    << " for cell " << cellI
                    << abort(FatalError);
            }
        }

        // Write split points into a point set
        pointSet pSet
        (
            mesh_,
            "splitPoints",
            labelHashSet(splitPointsToUnrefine_)
        );
        pSet.write();

        // Write split point cells into a cell set
        cellSet cSet
        (
            mesh_,
            "splitPointCells",
            splitPointsToUnrefine_.size()
        );

        forAll(splitPointsToUnrefine_, i)
        {
            // Get point cells and insert them into cell set
            const labelList& pCells = meshPointCells[splitPointsToUnrefine_[i]];

            forAll(pCells, j)
            {
                cSet.insert(pCells[j]);
            }
        }
        cSet.write();

        Pout<< "prismatic2DRefinement::setPrismatic2DUnrefinement"
            << "(polyTopoChange& ref)"
            << nl
            << "Writing " << pSet.size()
            << " points and "
            << cSet.size() << " cells for unrefinement to" << nl
            << "pointSet " << pSet.objectPath() << nl
            << "cellSet " << cSet.objectPath()
            << endl;
    }

    // Update refinementLevelIndicator for all cells that will be unrefined
    forAll(splitPointsToUnrefine_, i)
    {
        // Get point cells and mark them for unrefinement
        const labelList& pCells = meshPointCells[splitPointsToUnrefine_[i]];

        forAll(pCells, j)
        {
            refinementLevelIndicator_[pCells[j]] = UNREFINED;
        }
    }

    // Create lists needed by face remover
    labelList cellRegion;
    labelList cellRegionMaster;
    labelList facesToRemove;

    // Memory management
    {
        // Mark faces on empty patches to exclude them
        boolList faceOnEmptyPatch(mesh_.nFaces(), false);

        // Get boundary
        const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();

        // Loop through all patches
        forAll (boundaryMesh, patchI)
        {
            // Get current patch
            const polyPatch& curPatch = boundaryMesh[patchI];

            // Check whether this patch is emptyPolyPatch
            if (isA<emptyPolyPatch>(curPatch))
            {
                // Get start and end face labels
                const label startFaceI = curPatch.start();
                const label endFaceI = startFaceI + curPatch.size();

                // Mark all the faces and edges on the patch
                for (label faceI = startFaceI; faceI < endFaceI; ++faceI)
                {
                    // Mark face
                    faceOnEmptyPatch[faceI] = true;
                }
            }
        }


        // Collect split faces in the hash set, guess size to prevent excessive
        // resizing
        labelHashSet splitFaces(12*splitPointsToUnrefine_.size());

        // Get point faces
        const labelListList& meshPointFaces = mesh_.pointFaces();

        forAll(splitPointsToUnrefine_, i)
        {
            // Loop through all faces of this point and insert face index
            const labelList& pFaces = meshPointFaces[splitPointsToUnrefine_[i]];

            forAll(pFaces, j)
            {
                // Get face index
                const label& faceI = pFaces[j];

                if (!faceOnEmptyPatch[faceI])
                {
                    // Face is not on empty patch, insert it into hash set
                    splitFaces.insert(faceI);
                }
            }
        }

        // Check with faceRemover what faces will get removed. Note that this
        // can be more (but never less) than splitFaces provided.
        faceRemover_.compatibleRemoves
        (
            splitFaces.toc(),   // Pierced faces

            cellRegion,         // Region merged into (-1 for no region)
            cellRegionMaster,   // Master cell for region
            facesToRemove       // List of faces to be removed
        );

        if (facesToRemove.size() != splitFaces.size())
        {
            FatalErrorIn
            (
                "prismatic2DRefinement::setPrismatic2DUnrefinement"
                "(polyTopoChange& ref)"
            )   << "Either the initial set of split points to unrefine does not"
                << " seem to be consistent or there are no mid points of"
                << " refined cells."
                << abort(FatalError);
        }
    }

    // Insert all commands to combine cells
    faceRemover_.setRefinement
    (
        facesToRemove,
        cellRegion,
        cellRegionMaster,
        ref
    );
}


Foam::label Foam::prismatic2DRefinement::addFace
(
    polyTopoChange& ref,
    const label faceI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    // Set face information
    label patchID, zoneID, zoneFlip;
    meshTools::setFaceInfo(mesh_, faceI, patchID, zoneID, zoneFlip);

    // Set new face index to -1
    label newFaceI = -1;

    if ((nei == -1) || (own < nei))
    {
        // Ordering is ok, add the face
        newFaceI = ref.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                faceI,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
    else
    {
        // Ordering is flipped, reverse face and flip owner/neighbour
        newFaceI = ref.setAction
        (
            polyAddFace
            (
                newFace.reverseFace(),      // face
                nei,                        // owner
                own,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                faceI,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }

    return newFaceI;
}


Foam::label Foam::prismatic2DRefinement::addInternalFace
(
    polyTopoChange& ref,
    const label meshFaceI,
    const label meshPointI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    // Check whether this is an internal face
    if (mesh_.isInternalFace(meshFaceI))
    {
        return ref.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                meshFaceI,                  // master face for addition
                false,                      // flux flip
                -1,                         // patch for face
                -1,                         // zone for face
                false                       // face zone flip
            )
        );
    }
    else
    {
        // This is not an internal face. Add face out of nothing
        return ref.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                0,                          // master face for addition
                false,                      // flux flip
                -1,                         // patch for face
                -1,                         // zone for face
                false                       // face zone flip
            )
        );
    }
}


void Foam::prismatic2DRefinement::modifyFace
(
    polyTopoChange& ref,
    const label faceI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    // Set face inforomation
    label patchID, zoneID, zoneFlip;
    meshTools::setFaceInfo(mesh_, faceI, patchID, zoneID, zoneFlip);

    // Get owner/neighbour addressing and mesh faces
    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    const faceList& meshFaces = mesh_.faces();

    if
    (
        (own != owner[faceI])
     || (
            mesh_.isInternalFace(faceI)
         && (nei != neighbour[faceI])
        )
     || (newFace != meshFaces[faceI])
    )
    {
        // Either:
        // 1. Owner index does not correspond to mesh owner,
        // 2. Neighbour index does not correspond to mesh neighbour,
        // 3. New face does not correspond to mesh face
        // So we need to modify this face
        if ((nei == -1) || (own < nei))
        {
            // Ordering is ok, add the face
            ref.setAction
            (
                polyModifyFace
                (
                    newFace,            // modified face
                    faceI,              // label of face being modified
                    own,                // owner
                    nei,                // neighbour
                    false,              // face flip
                    patchID,            // patch for face
                    false,              // remove from zone
                    zoneID,             // zone for face
                    zoneFlip            // face flip in zone
                )
            );
        }
        else
        {
            // Ordering is flipped, reverse face and flip owner/neighbour
            ref.setAction
            (
                polyModifyFace
                (
                    newFace.reverseFace(),  // modified face
                    faceI,                  // label of face being modified
                    nei,                    // owner
                    own,                    // neighbour
                    false,                  // face flip
                    patchID,                // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }
    }
}


void Foam::prismatic2DRefinement::createInternalFaces
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& faceMidPoint,
    const labelList& faceAnchorLevel,
    const labelList& edgeMidPoint,
    const label cellI,

    polyTopoChange& ref
) const
{
    // Find in every face the cellLevel + 1 points (from edge subdivision)
    // and the anchor points

    // Get current cell and its level
    const cell& curCell = mesh_.cells()[cellI];
    const label& cLevel = cellLevel_[cellI];

    // Get mesh faces and face edges
    const faceList& meshFaces = mesh_.faces();
    const labelListList& meshFaceEdges = mesh_.faceEdges();

    // Get mesh cell points
    const labelListList& meshCellPoints = mesh_.cellPoints();

    // Map from edge mid to anchor points
    Map<edge> midPointToAnchors(24);
    // Map from edge mid to face mids
    Map<edge> midPointToFaceMids(24);

    // Running count of number of internal faces added so far
    label nFacesAdded = 0;

    // Loop through faces of the cell
    forAll(curCell, i)
    {
        // Get face index
        const label& faceI = curCell[i];

        // Get current face and its edges
        const face& f = meshFaces[faceI];
        const labelList& fEdges = meshFaceEdges[faceI];

        // We are on the cellI side of face f. The face will have 1 or n
        // cLevel points (where n is the number of points/edges of a face)
        // and lots of higher numbered ones

        // Index of face mid point
        label faceMidPointI = -1;

        // Get number of anchors for the face
        const label nAnchors = countAnchors(f, cLevel);

        if (nAnchors == 1)
        {
            // Only one anchor point. So the other side of the face has already
            // been split using cLevel + 1 and cLevel + 2 points

            // Find the one anchor
            label anchorFp = -1;

            // Loop through face points
            forAll(f, fp)
            {
                if (pointLevel_[f[fp]] <= cLevel)
                {
                    // Point level is smaller than cLevel + 1, this is the
                    // anchor point
                    anchorFp = fp;
                    break;
                }
            }

            // Now the face mid point is the second cLevel + 1 point
            label edgeMid = findLevel(f, f.fcIndex(anchorFp), true, cLevel + 1);
            label faceMid = findLevel(f, f.fcIndex(edgeMid), true, cLevel + 1);

            // Set face mid point index
            faceMidPointI = f[faceMid];
        }
        else
        {
            // There is no face middle yet but the face will be split. Set face
            // mid point index
            faceMidPointI = faceMidPoint[faceI];
        }


        // Now loop over all the anchors (might be just one) and store
        // the edge mids connected to it. storeMidPointInfo will collect
        // all the info and combine it all
        forAll(f, fp0)
        {
            // Get point index
            const label& point0 = f[fp0];

            if (pointLevel_[point0] <= cLevel)
            {
                // This is anchor point

                // Walk forward to cLevel + 1 or edgeMidPoint of this level
                label edgeMidPointI = -1;

                const label fp1 = f.fcIndex(fp0);

                if (pointLevel_[f[fp1]] <= cLevel)
                {
                    // Another anchor: edge will be split
                    const label& edgeI = fEdges[fp0];

                    edgeMidPointI = edgeMidPoint[edgeI];

                    // Sanity check
                    if (edgeMidPointI == -1)
                    {
                        const labelList& cPoints = meshCellPoints[cellI];

                        FatalErrorIn
                        (
                            "prismatic2DRefinement::createInternalFaces(...)"
                        )   << "cell:" << cellI << " cLevel:" << cLevel
                            << " cell points:" << cPoints
                            << " pointLevel:"
                            << IndirectList<label>(pointLevel_, cPoints)()
                            << " face:" << faceI
                            << " f:" << f
                            << " pointLevel:"
                            << IndirectList<label>(pointLevel_, f)()
                            << " faceAnchorLevel:" << faceAnchorLevel[faceI]
                            << " faceMidPoint:" << faceMidPoint[faceI]
                            << " faceMidPointI:" << faceMidPointI
                            << " fp:" << fp0
                            << abort(FatalError);
                    }
                }
                else
                {
                    // Search forward in face to clevel + 1
                    const label edgeMid = findLevel(f, fp1, true, cLevel + 1);

                    edgeMidPointI = f[edgeMid];
                }

                label newFaceI = storeMidPointInfo
                (
                    cellAnchorPoints,
                    cellAddedCells,
                    cellMidPoint,
                    edgeMidPoint,

                    cellI,
                    faceI,
                    true,                   // mid point after anchor
                    edgeMidPointI,          // edgemid
                    point0,                 // anchor
                    faceMidPointI,

                    midPointToAnchors,
                    midPointToFaceMids,
                    ref
                );

                if (newFaceI != -1)
                {
                    ++nFacesAdded;
                }


                // Now walk backward

                label fpMin1 = f.rcIndex(fp0);

                if (pointLevel_[f[fpMin1]] <= cLevel)
                {
                    // Another anchor: edge will be split
                    const label& edgeI = fEdges[fpMin1];

                    edgeMidPointI = edgeMidPoint[edgeI];

                    // Sanity check
                    if (edgeMidPointI == -1)
                    {
                        const labelList& cPoints = meshCellPoints[cellI];

                        FatalErrorIn
                        (
                            "prismatic2DRefinement::createInternalFaces(...)"
                        )
                            << "cell:" << cellI << " cLevel:" << cLevel
                            << " cell points:" << cPoints
                            << " pointLevel:"
                            << IndirectList<label>(pointLevel_, cPoints)()
                            << " face:" << faceI
                            << " f:" << f
                            << " pointLevel:"
                            << IndirectList<label>(pointLevel_, f)()
                            << " faceAnchorLevel:" << faceAnchorLevel[faceI]
                            << " faceMidPoint:" << faceMidPoint[faceI]
                            << " faceMidPointI:" << faceMidPointI
                            << " fp:" << fp0
                            << abort(FatalError);
                    }
                }
                else
                {
                    // Search back in face to clevel + 1
                    const label edgeMid = findLevel(f, fpMin1, false, cLevel + 1);

                    edgeMidPointI = f[edgeMid];
                }

                newFaceI = storeMidPointInfo
                (
                    cellAnchorPoints,
                    cellAddedCells,
                    cellMidPoint,
                    edgeMidPoint,

                    cellI,
                    faceI,
                    false,                  // mid point before anchor
                    edgeMidPointI,          // edgemid
                    point0,                 // anchor
                    faceMidPointI,

                    midPointToAnchors,
                    midPointToFaceMids,
                    ref
                );

                if (newFaceI != -1)
                {
                    ++nFacesAdded;
                }
            } // End for this anchor point
        } // End for all face points
    } // End for all cell faces
}


Foam::label Foam::prismatic2DRefinement::getAnchorCell
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label cellI,
    const label faceI,
    const label pointI
) const
{
    // Check whether pointI is an anchor of cellI. If it is not, check whether
    // any other point on the face is an anchor for the cell

    if (cellAnchorPoints[cellI].size() > 0)
    {
        const label index = findIndex(cellAnchorPoints[cellI], pointI);

        if (index != -1)
        {
            return cellAddedCells[cellI][index];
        }


        // pointI is not an anchor for the cell. Maybe we already refined the
        // face so check all the face vertices
        const face& f = mesh_.faces()[faceI];

        forAll(f, fp)
        {
            const label index = findIndex(cellAnchorPoints[cellI], f[fp]);

            if (index != -1)
            {
                return cellAddedCells[cellI][index];
            }
        }

        // Problem: the point does not seem to be an anchor for cell

        // Pick up points of the cell
        const labelList& cPoints = mesh_.cellPoints()[cellI];

        Perr<< "cell: " << cellI << ", points: " << endl;
        forAll(cPoints, i)
        {
            const label pointI = cPoints[i];

            Perr<< "    " << pointI << " coord: " << mesh_.points()[pointI]
                << nl;
        }

        Perr<< "cell: " << cellI << " anchorPoints: " << cellAnchorPoints[cellI]
            << endl;

        FatalErrorIn("prismatic2DRefinement::getAnchorCell(...)")
            << "Could not find point " << pointI
            << " in the anchorPoints for cell " << cellI << endl
            << "Does your original mesh obey the 2:1 constraint and"
            << " did you use consistentRefinement to make your cells to refine"
            << " obey this constraint as well?"
            << abort(FatalError);

        return -1;
    }
    else
    {
        return cellI;
    }
}


void Foam::prismatic2DRefinement::setNewFaceNeighbours
(
    const HashTable
    <
        label,
        Pair<label>,
        Hash<FixedList<label, 2> >
    >& pointCellToAddedCellMap,
    const labelListList& cellAddedCells,
    const label& faceI,
    const label& pointI,

    label& own,
    label& nei
) const
{
    // Get anchor cell for this anchor point on owner side
    const label ownCellI = mesh_.faceOwner()[faceI];

    // Get cell added cells for owner cell
    const labelList& cAddedOwn = cellAddedCells[ownCellI];

    // If the cell added cells list is not empty, return the necessary cell
    // index fetched with pointCellToAddedCellMap
    const Pair<label> pointOwnerCellPair(pointI, ownCellI);

    if (cAddedOwn.empty())
    {
        // Cell added cells is empty for owner, fetch original owner
        own = ownCellI;
    }
    else if (!pointCellToAddedCellMap.found(pointOwnerCellPair))
    {
        // Point-cell pair not found, meaning that we need to search for the
        // closest point which has lower level than this point. This happens if
        // we are splitting a face that has already been split but has different
        // owner/neighbour levels (e.g. owner was refined, but the neighbour was
        // not in the previous time step). It is possible that we end up with
        // splitting this face with point levels e.g. (1 1 0 0). Therefore, we
        // need to search for the point with minimum level on the edges sharing
        // this point.

        // Find point in the local faces addressing
        const face& f = mesh_.faces()[faceI];
        const label fpI = findIndex(f, pointI);

        if (fpI == -1)
        {
            FatalErrorIn("void prismatic2DRefinement::setNewFaceNeighbours(...)")
                << "Point: " << pointI << " not found in face: " << f
                << ", with face index: " << faceI
                << nl
                << "Point must belong to the face."
                << nl
                << "Error encounted when dealing with owner cell: "
                << ownCellI
                << abort(FatalError);
        }

        // Find the global point index with minimum edge connected level
        const label anchorPointI = findMinEdgeConnectedLevel
        (
            fpI,                      // current face point
            faceI,                    // face index
            f,                        // face
            mesh_.faceEdges()[faceI], // face edges
            mesh_.edges()             // mesh edges
        );

        // If the point is the same as pointI, we did not find any valid point
        if (anchorPointI == pointI)
        {
            FatalErrorIn("void prismatic2DRefinement::setNewFaceNeighbours(...)")
                << "Could not find different adjacent anchor point."
                << nl
                << "pointI: " << pointI << " faceI: " << faceI
                << "Error encounted when dealing with owner cell: "
                << ownCellI
                << abort(FatalError);
        }

        // Set owner cell
        own = cAddedOwn
        [
            pointCellToAddedCellMap[Pair<label>(anchorPointI, ownCellI)]
        ];
    }
    else
    {
        // Cell added cells is not empty and the mapping is found. Set owner
        // from mapping
        own = cAddedOwn[pointCellToAddedCellMap[pointOwnerCellPair]];
    }

    if (mesh_.isInternalFace(faceI))
    {
        // Get anchor cell for this anchor point on neighbour side
        const label neiCellI = mesh_.faceNeighbour()[faceI];

        // Get cell added cells for neighbour cell
        const labelList& cAddedNei = cellAddedCells[neiCellI];

        // If the cell added cells list is not empty, return the necessary cell
        // index fetched with pointCellToAddedCellMap

        const Pair<label> pointNeighbourCellPair(pointI, neiCellI);

        if (cAddedNei.empty())
        {
            // Cell added cells is empty for neighbour, fetch original neighbour
            nei = neiCellI;
        }
        else if (!pointCellToAddedCellMap.found(pointNeighbourCellPair))
        {
            // Point-cell pair not found, meaning that we need to search for the
            // closest point which has lower level than this point. This happens
            // if we are splitting a face that has already been split but has
            // different owner/neighbour levels (e.g. owner was refined, but the
            // neighbour was not in the previous time step). It is possible that
            // we end up with splitting this face with point levels e.g. (1 1 0
            // 0). Therefore, we need to search for the point with minimum level
            // on the edges sharing this point.

            // Find point in the local faces addressing
            const face& f = mesh_.faces()[faceI];
            const label fpI = findIndex(f, pointI);

            if (fpI == -1)
            {
                FatalErrorIn("void prismatic2DRefinement::setNewFaceNeighbours(...)")
                    << "Point: " << pointI << " not found in face: " << f
                    << ", with face index: " << faceI
                    << nl
                    << "Point must belong to the face."
                    << nl
                    << "Error encounted when dealing with neighbour cell: "
                    << neiCellI
                    << abort(FatalError);
            }

            // Find the global point index with minimum edge connected level
            const label anchorPointI = findMinEdgeConnectedLevel
            (
                fpI,                      // current face point
                faceI,                    // face index
                f,                        // face
                mesh_.faceEdges()[faceI], // face edges
                mesh_.edges()             // mesh edges
            );

            // If the point is the same as pointI, we did not find any valid point
            if (anchorPointI == pointI)
            {
                FatalErrorIn("void prismatic2DRefinement::setNewFaceNeighbours(...)")
                    << "Could not find different adjacent anchor point."
                    << nl
                    << "pointI: " << pointI << " faceI: " << faceI
                    << "Error encounted when dealing with neighbour cell: "
                    << neiCellI
                    << abort(FatalError);
            }

            // Set owner cell
            nei = cAddedNei
            [
                pointCellToAddedCellMap[Pair<label>(anchorPointI, neiCellI)]
            ];
        }
        else
        {
            // Cell added cells is not empty and the mapping is found. Set
            // neighbour from mapping
            nei = cAddedNei[pointCellToAddedCellMap[pointNeighbourCellPair]];
        }
    }
    else
    {
        // Boundary face: set neighbour to -1
        nei = -1;
    }
}


void Foam::prismatic2DRefinement::walkFaceToMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label faceI,
    const label startFp,
    dynamicLabelList& faceVerts
) const
{
    // Get the face and its edges
    const face& f = mesh_.faces()[faceI];
    const labelList& fEdges = mesh_.faceEdges()[faceI];

    label fp = startFp;

    // Starting from fp store all (1 or 2) vertices until where the face
    // gets split
    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] > -1)
        {
            // Edge is split, append its mid point
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (pointLevel_[f[fp]] <= cLevel)
        {
            // Found next anchor. Already appended the split point just above
            return;
        }
        else if (pointLevel_[f[fp]] == cLevel + 1)
        {
            // Mid level, append and return
            faceVerts.append(f[fp]);

            return;
        }
        else if (pointLevel_[f[fp]] == cLevel + 2)
        {
            // Store and continue to cLevel + 1
            faceVerts.append(f[fp]);
        }
    }
}


void Foam::prismatic2DRefinement::walkFaceFromMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label faceI,
    const label startFp,
    dynamicLabelList& faceVerts
) const
{
    // Get the face and its edges
    const face& f = mesh_.faces()[faceI];
    const labelList& fEdges = mesh_.faceEdges()[faceI];

    label fp = f.rcIndex(startFp);

    while (true)
    {
        if (pointLevel_[f[fp]] <= cLevel)
        {
            // Anchor point, break
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel + 1)
        {
            // Mid level, append and break
            faceVerts.append(f[fp]);
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel + 2)
        {
            // Continue to cLevel + 1
        }
        fp = f.rcIndex(fp);
    }

    // Store
    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] > -1)
        {
            // Edge is split, append its mid point
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (fp == startFp)
        {
            break;
        }
        faceVerts.append(f[fp]);
    }
}


Foam::label Foam::prismatic2DRefinement::findMinLevel(const labelList& f) const
{
    // Initialise minimum level to large value
    label minLevel = labelMax;

    // Initialise point label at which min level is reached to -1
    label pointIMin = -1;

    forAll(f, fp)
    {
        const label& level = pointLevel_[f[fp]];

        if (level < minLevel)
        {
            minLevel = level;
            pointIMin = fp;
        }
    }

    return pointIMin;
}


Foam::label Foam::prismatic2DRefinement::findMinEdgeConnectedLevel
(
    const label& fpI,
    const label& faceI,
    const face& f,
    const labelList& fEdges,
    const edgeList& meshEdges
) const
{
    // Get point index and initialize anchor point
    const label& pointI = f[fpI];

    label anchorPointI = pointI;

    // Check the other point on edge starting with pointI
    const label& edgeIndexAfterPoint = fEdges[fpI];
    const edge& edgeAfter = meshEdges[edgeIndexAfterPoint];

    // Get other point on edge after
    const label& pointAfter = edgeAfter.otherVertex(pointI);

    if (pointLevel_[pointAfter] < pointLevel_[anchorPointI])
    {
        anchorPointI = pointAfter;
    }

    // Check the other point on edge ending with pointI
    const label& edgeIndexBeforePoint = fEdges[f.rcIndex(fpI)];
    const edge& edgeBefore = meshEdges[edgeIndexBeforePoint];

    // Get other point on edge before
    const label& pointBefore = edgeBefore.otherVertex(pointI);

    if (pointLevel_[pointBefore] < pointLevel_[anchorPointI])
    {
        anchorPointI = pointBefore;
    }

    return anchorPointI;
}


Foam::label Foam::prismatic2DRefinement::findMaxLevel(const labelList& f) const
{
    // Initialise  maximum level to small value
    label maxLevel = labelMin;

    // Initialise point label at which max level is reached to -1
    label pointIMax = -1;

    forAll(f, fp)
    {
        const label& level = pointLevel_[f[fp]];

        if (level > maxLevel)
        {
            maxLevel = level;
            pointIMax = fp;
        }
    }

    return pointIMax;
}


Foam::label Foam::prismatic2DRefinement::countAnchors
(
    const labelList& f,
    const label anchorLevel
) const
{
    label nAnchors = 0;

    forAll(f, fp)
    {
        if (pointLevel_[f[fp]] <= anchorLevel)
        {
            ++nAnchors;
        }
    }
    return nAnchors;
}


void Foam::prismatic2DRefinement::addFaceMids
(
    const labelList& faceMidPoint,
    const boolList& faceOnEmptyPatch,
    const label& faceI,
    const label& cellI,
    face& newFace
) const
{
    // b) faceMidPoint for this face
    newFace[1] = faceMidPoint[faceI];

    // c) faceMidPoint for the face on the other side
    // The other face is uniquely defined as the other face of the same cell
    // which is on empty patch

    // Get the cell
    const cell& cFaces = mesh_.cells()[cellI];

    // Loop through cell faces
    forAll(cFaces, i)
    {
        // Get face index
        const label& faceJ = cFaces[i];

        if (faceOnEmptyPatch[faceJ] && (faceI != faceJ))
        {
            // This is the face we're looking for, add its
            // midpoint and double check if it is valid
            if (faceMidPoint[faceJ] > -1)
            {
                newFace[2] = faceMidPoint[faceJ];
                break;
            }
            else
            {
                FatalErrorIn
                (
                    "void prismatic2DRefinement::addFaceMids(...)"
                )   << "Other face: "
                    << faceJ
                    << " has not been selected for splitting,"
                    << " while the face on original side: "
                    << faceI
                    <<" has been selected."
                    << abort(FatalError);
            }
        } // End if this is our face
    } // End for all cell faces
}


Foam::label Foam::prismatic2DRefinement::findLevel
(
    const face& f,
    const label startFp,
    const bool searchForward,
    const label wantedLevel
) const
{
    label fp = startFp;

    forAll(f, i)
    {
        label pointI = f[fp];

        if (pointLevel_[pointI] < wantedLevel)
        {
            FatalErrorIn("prismatic2DRefinement::findLevel(...)")
                << "face:" << f
                << " level:" << IndirectList<label>(pointLevel_, f)()
                << " startFp:" << startFp
                << " wantedLevel:" << wantedLevel
                << abort(FatalError);
        }
        else if (pointLevel_[pointI] == wantedLevel)
        {
            return fp;
        }

        if (searchForward)
        {
            fp = f.fcIndex(fp);
        }
        else
        {
            fp = f.rcIndex(fp);
        }
    }

    FatalErrorIn("prismatic2DRefinement::findLevel(...)")
        << "face:" << f
        << " level:" << IndirectList<label>(pointLevel_, f)()
        << " startFp:" << startFp
        << " wantedLevel:" << wantedLevel
        << abort(FatalError);

    return -1;
}


Foam::label Foam::prismatic2DRefinement::storeMidPointInfo
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& edgeMidPoint,
    const label cellI,
    const label faceI,
    const bool faceOrder,
    const label edgeMidPointI,
    const label anchorPointI,
    const label faceMidPointI,

    Map<edge>& midPointToAnchors,
    Map<edge>& midPointToFaceMids,
    polyTopoChange& ref
) const
{
    // A single internal face is added per edge inbetween anchor points,
    // i.e. one face per midPoint between anchor points. The information is
    // stored on the midPoint and if we have enough information (finished
    // collecting two anchors and two face mid points), we add the face.
    // Note that this member function can get called anywhere from
    // two times (two unrefined faces) to four times (two refined faces) so
    // the first call that adds the information creates the face


    // See if need to store anchors
    bool changed = false;
    bool haveTwoAnchors = false;

    Map<edge>::iterator edgeMidFnd = midPointToAnchors.find(edgeMidPointI);

    if (edgeMidFnd == midPointToAnchors.end())
    {
        midPointToAnchors.insert(edgeMidPointI, edge(anchorPointI, -1));
    }
    else
    {
        edge& e = edgeMidFnd();

        if (anchorPointI != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = anchorPointI;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoAnchors = true;
        }
    }

    bool haveTwoFaceMids = false;

    Map<edge>::iterator faceMidFnd = midPointToFaceMids.find(edgeMidPointI);

    if (faceMidFnd == midPointToFaceMids.end())
    {
        midPointToFaceMids.insert(edgeMidPointI, edge(faceMidPointI, -1));
    }
    else
    {
        edge& e = faceMidFnd();

        if (faceMidPointI != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = faceMidPointI;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoFaceMids = true;
        }
    }

    // Check if this call of storeMidPointInfo is the one that completed all
    // the nessecary information

    if (changed && haveTwoAnchors && haveTwoFaceMids)
    {
        const edge& anchors = midPointToAnchors[edgeMidPointI];
        const edge& faceMids = midPointToFaceMids[edgeMidPointI];

        label otherFaceMidPointI = faceMids.otherVertex(faceMidPointI);

        // Create face consistent with anchorI being the owner.
        // Note that the edges between the edge mid point and the face mids
        // might be marked for splitting. Note that these edge splits cannot
        // be between cell mid and face mids

        dynamicLabelList newFaceVerts(4);
        if (faceOrder == (mesh_.faceOwner()[faceI] == cellI))
        {
            newFaceVerts.append(faceMidPointI);

            // Check and insert edge split if any
            insertEdgeSplit
            (
                edgeMidPoint,
                faceMidPointI,  // edge between faceMid and
                edgeMidPointI,  // edgeMid
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointI,
                otherFaceMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(otherFaceMidPointI);
            newFaceVerts.append(cellMidPoint[cellI]);
        }
        else
        {
            newFaceVerts.append(otherFaceMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                otherFaceMidPointI,
                edgeMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointI,
                faceMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(faceMidPointI);
            newFaceVerts.append(cellMidPoint[cellI]);
        }

        face newFace;
        newFace.transfer(newFaceVerts.shrink());
        newFaceVerts.clear();

        label anchorCell0 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            cellI,
            faceI,
            anchorPointI
        );
        label anchorCell1 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            cellI,
            faceI,
            anchors.otherVertex(anchorPointI)
        );

        // Get mesh points
        const pointField& meshPoints = mesh_.points();

        label own, nei;
        point ownPt, neiPt;

        if (anchorCell0 < anchorCell1)
        {
            own = anchorCell0;
            nei = anchorCell1;

            ownPt = meshPoints[anchorPointI];
            neiPt = meshPoints[anchors.otherVertex(anchorPointI)];

        }
        else
        {
            own = anchorCell1;
            nei = anchorCell0;
            newFace = newFace.reverseFace();

            ownPt = meshPoints[anchors.otherVertex(anchorPointI)];
            neiPt = meshPoints[anchorPointI];
        }

        if (debug)
        {
            point ownPt, neiPt;

            if (anchorCell0 < anchorCell1)
            {
                ownPt = meshPoints[anchorPointI];
                neiPt = meshPoints[anchors.otherVertex(anchorPointI)];
            }
            else
            {
                ownPt = meshPoints[anchors.otherVertex(anchorPointI)];
                neiPt = meshPoints[anchorPointI];
            }

            checkInternalOrientation
            (
                ref,
                cellI,
                faceI,
                ownPt,
                neiPt,
                newFace
            );
        }

        return addInternalFace
        (
            ref,
            faceI,
            anchorPointI,
            newFace,
            own,
            nei
        );
    }
    else
    {
        return -1;
    }
}


void Foam::prismatic2DRefinement::insertEdgeSplit
(
    const labelList& edgeMidPoint,
    const label p0,
    const label p1,
    dynamicLabelList& verts
) const
{
    // Get number of points
    const label nPoints = mesh_.nPoints();

    if (p0 < nPoints && p1 < nPoints)
    {
        label edgeI = meshTools::findEdge(mesh_, p0, p1);

        if (edgeI != -1 && edgeMidPoint[edgeI] != -1)
        {
            verts.append(edgeMidPoint[edgeI]);
        }
    }
}


void Foam::prismatic2DRefinement::checkNewFaceOrientation
(
    polyTopoChange& ref,
    const label& faceI,
    const face& newFace
) const
{
    // Get mesh cell centres
    const vectorField& meshCellCentres = mesh_.cellCentres();

    if (mesh_.isInternalFace(faceI))
    {
        // Get old owner/neighbour indices
        const label oldOwn = mesh_.faceOwner()[faceI];
        const label oldNei = mesh_.faceNeighbour()[faceI];

        // Print info only with deep debug level
        if (debug > 1)
        {
            Pout<< "Split infternal face: " << faceI
                << ", into quad: " << newFace << nl
                << "owner: " << oldOwn
                << ", neighbour: " << oldNei << endl;
        }

        checkInternalOrientation
        (
            ref,
            oldOwn,
            faceI,
            meshCellCentres[oldOwn],
            meshCellCentres[oldNei],
            newFace
        );
    }
    else
    {
        // Get face centres and old owner
        const vectorField& meshFaceCentres = mesh_.faceCentres();
        const label oldOwn = mesh_.faceOwner()[faceI];

        // Print info only with deep debug level
        if (debug > 1)
        {
            Pout<< "Split boundary face: " << faceI
                << ", into quad: " << newFace << nl
                << "owner: " << oldOwn << endl;
        }

        checkBoundaryOrientation
        (
            ref,
            oldOwn,
            faceI,
            meshCellCentres[oldOwn],
            meshFaceCentres[faceI],
            newFace
        );
    }
}


void Foam::prismatic2DRefinement::checkInternalOrientation
(
    polyTopoChange& ref,
    const label cellI,
    const label faceI,
    const point& ownPt,
    const point& neiPt,
    const face& newFace
) const
{
    const face compactFace(identity(newFace.size()));

    // Get list of polyAddPoint objects
    const DynamicList<polyAddPoint>& polyAddedPoints(ref.addedPoints());

    // Create a list of all points
    const label nOrigPoints = mesh_.nPoints();
    pointField allPoints(nOrigPoints + polyAddedPoints.size());

    // Set ordinary points first
    const pointField& meshPoints = mesh_.points();
    forAll (meshPoints, i)
    {
        allPoints[i] = meshPoints[i];
    }

    // Set newly added points next
    forAll (polyAddedPoints, i)
    {
        allPoints[i + nOrigPoints] = polyAddedPoints[i].newPoint();
    }

    // Get compact points
    const pointField compactPoints
    (
        IndirectList<point>(allPoints, newFace)()
    );

    const vector n(compactFace.normal(compactPoints));
    const vector dir(neiPt - ownPt);

    // Check orientation error
    if ((dir & n) < 0)
    {
        FatalErrorIn
        (
            "prismatic2DRefinement::checkInternalOrientation(...)"
        )
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace << endl
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " neiPt:" << neiPt
            << abort(FatalError);
    }

    // Note: report significant non-orthogonality error
    const scalar severeNonOrthogonalityThreshold =
      ::cos
        (
            primitiveMesh::nonOrthThreshold_()/180.0*mathematicalConstant::pi
        );

    const vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    // Note: normal vector already normalised
    const scalar dDotN = (fcToOwn & n)/(mag(fcToOwn) + VSMALL);

    if (dDotN > severeNonOrthogonalityThreshold)
    {
        WarningIn
        (
            "prismatic2DRefinement::checkInternalOrientation(...)"
        )
            << "Detected severely non-orthogonal face with non-orthogonality: "
            << ::acos(dDotN)/mathematicalConstant::pi*180.0
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace << endl
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " neiPt:" << neiPt
            << " dDotN:" << dDotN
            << endl;
    }
}


void Foam::prismatic2DRefinement::checkBoundaryOrientation
(
    polyTopoChange& ref,
    const label cellI,
    const label faceI,
    const point& ownPt,
    const point& boundaryPt,
    const face& newFace
) const
{
    const face compactFace(identity(newFace.size()));

    // Get list of polyAddPoint objects
    const DynamicList<polyAddPoint>& polyAddedPoints(ref.addedPoints());

    // Create a list of all points
    const label nOrigPoints = mesh_.nPoints();
    pointField allPoints(nOrigPoints + polyAddedPoints.size());

    // Set ordinary points first
    const pointField& meshPoints = mesh_.points();
    forAll (meshPoints, i)
    {
        allPoints[i] = meshPoints[i];
    }

    // Set newly added points next
    forAll (polyAddedPoints, i)
    {
        allPoints[i + nOrigPoints] = polyAddedPoints[i].newPoint();
    }

    // Get compact points
    const pointField compactPoints
    (
        IndirectList<point>(allPoints, newFace)()
    );

    const vector n(compactFace.normal(compactPoints));
    const vector dir(boundaryPt - ownPt);

    // Check orientation error
    if ((dir & n) < 0)
    {
        FatalErrorIn
        (
            "prismatic2DRefinement::checkBoundaryOrientation(...)"
        )   << "Detected invalid orientation of the boundary face." << nl
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << abort(FatalError);
    }

    // Note: report significant non-orthogonality error
    const scalar severeNonOrthogonalityThreshold =
      ::cos
        (
            primitiveMesh::nonOrthThreshold_()/180.0*mathematicalConstant::pi
        );

    const vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    // Note: normal vector already normalised
    const scalar dDotN = (fcToOwn & n)/(mag(fcToOwn) + VSMALL);

    if (dDotN > severeNonOrthogonalityThreshold)
    {
        WarningIn
        (
            "prismatic2DRefinement::checkBoundaryOrientation(...)"
        )
            << "Detected severely non-orthogonal face with non-orthogonality: "
            << ::acos(dDotN)/mathematicalConstant::pi*180.0
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << " dDotN:" << dDotN
            << endl;
    }
}


Foam::label Foam::prismatic2DRefinement::faceConsistentRefinement
(
    boolList& cellsToRefine
) const
{
    // Count number of cells that will be added
    label nAddCells = 0;

    // Get necessary mesh data
    const label nFaces = mesh_.nFaces();
    const label nInternalFaces = mesh_.nInternalFaces();

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    // Loop through internal faces and check consistency
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        // Get owner and neighbour labels
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        // Get owner and neighbour cell levels
        // Note: If the cell is marked for refinement, the level is current
        // level + 1, otherwise it is equal to the current level
        const label ownLevel =
            cellsToRefine[own] ? cellLevel_[own] + 1 : cellLevel_[own];
        const label neiLevel =
            cellsToRefine[nei] ? cellLevel_[nei] + 1 : cellLevel_[nei];

        if (ownLevel > (neiLevel + 1))
        {
            // Owner level is higher than neighbour level + 1, neighbour must be
            // marked for refinement
            cellsToRefine[nei] = true;
            ++nAddCells;
        }
        else if (neiLevel > (ownLevel + 1))
        {
            // Neighbour level is higher than owner level + 1, owner must be
            // marked for refinement
            cellsToRefine[own] = true;
            ++nAddCells;
        }
    }

    // Create owner level for boundary faces to prepare for swapping on coupled
    // boundaries
    labelList ownLevel(nFaces - nInternalFaces);
    forAll (ownLevel, i)
    {
        // Get owner of the face and update owner cell levels
        const label& own = owner[i + nInternalFaces];
        ownLevel[i] =
            cellsToRefine[own] ? cellLevel_[own] + 1 : cellLevel_[own];
    }

    // Swap boundary face lists (coupled boundary update)
    syncTools::swapBoundaryFaceList(mesh_, ownLevel, false);

    // Note: now the ownLevel list actually contains the neighbouring level
    // (from the other side), use alias (reference) for clarity from now on
    const labelList& neiLevel = ownLevel;

    // Loop through boundary faces
    forAll (neiLevel, i)
    {
        // Get owner of the face and owner level
        const label& own = owner[i + nInternalFaces];
        const label curOwnLevel =
            cellsToRefine[own] ? cellLevel_[own] + 1 : cellLevel_[own];

        // Note: we are using more stringent 1:1 consistency across coupled
        // boundaries in order to simplify handling of edge based consistency
        // checks for parallel runs
        if (neiLevel[i] > curOwnLevel)
        {
            // Neighbour level is higher than owner level, owner must be
            // marked for refinement
            cellsToRefine[own] = true;
            ++nAddCells;
        }

        // Note: other possibility (that owner level is higher than neighbour
        // level) is taken into account on the other side automatically
    }

    // Return number of added cells
    return nAddCells;
}


Foam::label Foam::prismatic2DRefinement::edgeConsistentRefinement
(
    boolList& cellsToRefine
) const
{
    // Count number of cells that will be added
    label nAddCells = 0;

    // Algorithm: loop over all edges and visit all unique cell pairs sharing
    // this particular edge. Then, ensure 2:1 edge consistency by marking
    // cell with lower level for refinement

    // Get edge cells
    const labelListList& meshEdgeCells = mesh_.edgeCells();

    // Loop through all mesh edges
    forAll (meshEdgeCells, edgeI)
    {
        // Get current edge cells
        const labelList& curEdgeCells = meshEdgeCells[edgeI];

        // Loop through all edge cells
        forAll (curEdgeCells, i)
        {
            // Get first cell index
            const label& cellI = curEdgeCells[i];

            // Loop through remaining edge cells
            for (label j = i + 1; j < curEdgeCells.size(); ++j)
            {
                // Get second cell index
                const label& cellJ = curEdgeCells[j];

                // Get levels of the two cells. If the cell is marked for
                // refinement, the level is current level + 1, otherwise it is
                // equal to the current level

                // Note: cellsToRefine flag for both cellI and cellJ might
                // change, this is why we need to recalculate cellI level here
                const label cellILevel =
                    cellsToRefine[cellI]
                  ? cellLevel_[cellI] + 1
                  : cellLevel_[cellI];

                const label cellJLevel =
                    cellsToRefine[cellJ]
                  ? cellLevel_[cellJ] + 1
                  : cellLevel_[cellJ];

                if (cellILevel > cellJLevel + 1)
                {
                    // Level of cellI is higher than level of cellJ + 1, cellJ
                    // must be marked for refinement
                    cellsToRefine[cellJ] = true;
                    ++nAddCells;
                }
                else if (cellJLevel > cellILevel + 1)
                {
                    // Level of cellJ is higher than level of cellI + 1, cellI
                    // must be marked for refinement
                    cellsToRefine[cellI] = true;
                    ++nAddCells;
                }
            }
        }
    }

    // Note: in order to avoid very difficult and time-consuming parallelisation
    // of edge cell connectivity and edge cell values, we enforce a more
    // stringent face-based consistency across processor boundaries. Basically,
    // if a face-based consistency of 1:1 (not 2:1 as for ordinary faces) is
    // ensured, the edge-based consistency becomes a local operation (I'm not
    // 100% sure to be honest since there are countless variants when dealing
    // with arbitrary prismatic2D cells).
    // See faceConsistentRefinement for details. VV, 17/Apr/2018.

    // Return number of added cells
    return nAddCells;
}


Foam::label Foam::prismatic2DRefinement::faceConsistentUnrefinement
(
    boolList& cellsToUnrefine
) const
{
    // Count number of removed cells from unrefinement
    label nRemCells = 0;

    // Get necessary mesh data
    const label nFaces = mesh_.nFaces();
    const label nInternalFaces = mesh_.nInternalFaces();

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    // Loop through internal faces and check consistency
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        // Get owner and neighbour labels
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        // Get owner and neighbour cell levels
        // Note: If the cell is marked for unrefinement, the level is current
        // level - 1, otherwise it is equal to the current level
        const label ownLevel =
            cellsToUnrefine[own] ? cellLevel_[own] - 1 : cellLevel_[own];
        const label neiLevel =
            cellsToUnrefine[nei] ? cellLevel_[nei] - 1 : cellLevel_[nei];

        if (ownLevel < (neiLevel - 1))
        {
            // Owner level is smaller than neighbour level - 1, we must not
            // unrefine owner

            // Check whether the cell has not been marked for unrefinement
            if (!cellsToUnrefine[own])
            {
                FatalErrorIn
                (
                    "label prismatic2DRefinement::faceConsistentUnrefinement"
                    "(boolList& cellsToUnrefine)"
                )   << "Cell not marked for unrefinement, indicating a"
                    << " previous unnoticed problem with unrefinement."
                    << nl
                    << "Owner: " << own << ", neighbour: " << nei
                    << nl
                    << "Owner level: " << ownLevel
                    << ", neighbour level: " << neiLevel << nl
                    << "This is probably because the refinement and "
                    << "unrefinement regions are very close." << nl
                    << "Try increasing nUnrefinementBufferLayers. "
                    << abort(FatalError);
            }

            cellsToUnrefine[own] = false;
            ++nRemCells;
        }
        else if (neiLevel < (ownLevel - 1))
        {
            // Neighbour level is smaller than owner level - 1, we must not
            // unrefine neighbour

            // Check whether the cell has not been marked for unrefinement
            if (!cellsToUnrefine[nei])
            {
                FatalErrorIn
                (
                    "label prismatic2DRefinement::faceConsistentUnrefinement"
                    "(boolList& cellsToUnrefine)"
                )   << "Cell not marked for unrefinement, indicating a"
                    << " previous unnoticed problem with unrefinement."
                    << nl
                    << "Owner: " << own << ", neighbour: " << nei
                    << nl
                    << "Owner level: " << ownLevel
                    << ", neighbour level: " << neiLevel << nl
                    << "This is probably because the refinement and "
                    << "unrefinement regions are very close." << nl
                    << "Try increasing nUnrefinementBufferLayers. "
                    << abort(FatalError);
            }

            cellsToUnrefine[nei] = false;
            ++nRemCells;
        }
    }

    // Create owner level for boundary faces to prepare for swapping on coupled
    // boundaries
    labelList ownLevel(nFaces - nInternalFaces);
    forAll (ownLevel, i)
    {
        // Get owner of the face and update owner cell levels
        const label& own = owner[i + nInternalFaces];
        ownLevel[i] =
            cellsToUnrefine[own] ? cellLevel_[own] - 1 : cellLevel_[own];
    }

    // Swap boundary face lists (coupled boundary update)
    syncTools::swapBoundaryFaceList(mesh_, ownLevel, false);

    // Note: now the ownLevel list actually contains the neighbouring level
    // (from the other side), use alias (reference) for clarity from now on
    const labelList& neiLevel = ownLevel;

    // Loop through boundary faces
    forAll (neiLevel, i)
    {
        // Get owner of the face and owner level
        const label& own = owner[i + nInternalFaces];
        const label curOwnLevel =
            cellsToUnrefine[own] ? cellLevel_[own] - 1 : cellLevel_[own];

        // Note: we are using more stringent 1:1 consistency across coupled
        // boundaries in order to simplify handling of edge based consistency
        // checkes for parallel runs
        if (curOwnLevel < neiLevel[i])
        {
            // Owner level is smaller than neighbour level, we must not
            // unrefine owner

            // Check whether the cell has not been marked for unrefinement
            if (!cellsToUnrefine[own])
            {
                FatalErrorIn
                (
                    "label prismatic2DRefinement::faceConsistentUnrefinement"
                    "(boolList& cellsToUnrefine)"
                )   << "Boundary cell not marked for unrefinement, indicating a"
                    << " previous unnoticed problem with unrefinement."
                    << nl
                    << "Owner: " << own
                    << nl
                    << "Owner level: " << curOwnLevel
                    << ", neighbour level: " << neiLevel[i] << nl
                    << "This is probably because the refinement and "
                    << "unrefinement regions are very close." << nl
                    << "Try increasing nUnrefinementBufferLayers. "
                    << abort(FatalError);
            }

            cellsToUnrefine[own] = false;
            ++nRemCells;
        }

        // Note: other possibility (that neighbour level is smaller than owner
        // level) is taken into account on the other side automatically
    }

    // Return number of local cells removed from unrefinement
    return nRemCells;
}


Foam::label Foam::prismatic2DRefinement::edgeConsistentUnrefinement
(
    boolList& cellsToUnrefine
) const
{
    // Count number of cells that will be removed
    label nRemCells = 0;

    // Algorithm: loop over all edges and visit all unique cell pairs sharing
    // this particular edge. Then, ensure 2:1 edge consistency by protecting the
    // cell with lower level from unrefinement

    // Get edge cells
    const labelListList& meshEdgeCells = mesh_.edgeCells();

    // Loop through all mesh edges
    forAll (meshEdgeCells, edgeI)
    {
        // Get current edge cells
        const labelList& curEdgeCells = meshEdgeCells[edgeI];

        // Loop through all edge cells
        forAll (curEdgeCells, i)
        {
            // Get first cell index
            const label& cellI = curEdgeCells[i];

            // Loop through remaining edge cells
            for (label j = i + 1; j < curEdgeCells.size(); ++j)
            {
                // Get second cell index
                const label& cellJ = curEdgeCells[j];

                // Get levels of the two cells. If the cell is marked for
                // unrefinement, the level is current level - 1, otherwise it is
                // equal to the current level

                // Note: cellsToUnrefine flag for both cellI and cellJ might
                // change, this is why we need to recalculate cellI level here
                const label cellILevel =
                    cellsToUnrefine[cellI]
                  ? cellLevel_[cellI] - 1
                  : cellLevel_[cellI];

                const label cellJLevel =
                    cellsToUnrefine[cellJ]
                  ? cellLevel_[cellJ] - 1
                  : cellLevel_[cellJ];

                if (cellILevel < cellJLevel - 1)
                {
                    // Level of cellI is smaller than level of cellJ - 1, cellI
                    // must be protected from unrefinement

                    // Check whether the cell has not been marked for
                    // unrefinement
                    if (!cellsToUnrefine[cellI])
                    {
                        FatalErrorIn
                        (
                            "label prismatic2DRefinement::"
                            "edgeConsistentUnrefinement"
                            "(boolList& cellsToUnrefine)"
                        )   << "Cell not marked for unrefinement, indicating a"
                            << " previous unnoticed problem with unrefinement."
                            << nl
                            << "cellI: " << cellI << ", cellJ: " << cellJ
                            << nl
                            << "Level of cellI: " << cellILevel
                            << ", level of cellJ: " << cellJLevel << nl
                            << "This is probably because the refinement and "
                            << "unrefinement regions are very close." << nl
                            << "Try increasing nUnrefinementBufferLayers. "
                            << abort(FatalError);
                    }

                    cellsToUnrefine[cellI] = false;
                    ++nRemCells;
                }
                else if (cellJLevel < cellILevel - 1)
                {
                    // Level of cellJ is smaller than level of cellI - 1, cellJ
                    // must be protected from unrefinement

                    // Check whether the cell has not been marked for
                    // unrefinement
                    if (!cellsToUnrefine[cellJ])
                    {
                        FatalErrorIn
                        (
                            "label prismatic2DRefinement::"
                            "edgeConsistentUnrefinement"
                            "(boolList& cellsToUnrefine)"
                        )   << "Cell not marked for unrefinement, indicating a"
                            << " previous unnoticed problem with unrefinement."
                            << nl
                            << "cellI: " << cellI << ", cellJ: " << cellJ
                            << nl
                            << "Level of cellI: " << cellILevel
                            << ", level of cellJ: " << cellJLevel << nl
                            << "This is probably because the refinement and "
                            << "unrefinement regions are very close." << nl
                            << "Try increasing nUnrefinementBufferLayers. "
                            << abort(FatalError);
                    }

                    cellsToUnrefine[cellJ] = false;
                    ++nRemCells;
                }
            }
        }
    }

    // Note: in order to avoid very difficult and time-consuming parallelisation
    // of edge cell connectivity and edge cell values, we enforce a more
    // stringent face-based consistency across processor boundaries. Basically,
    // if a face-based consistency of 1:1 (not 2:1 as for ordinary faces) is
    // ensured, the edge-based consistency becomes a local operation (I'm not
    // 100% sure to be honest whether this is true all the time since there are
    // countless variants when dealing with arbitrary prismatic2D cells).
    // See faceConsistentRefinement for details. VV, 3/Apr/2018.

    // Return number of removed cells
    return nRemCells;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::prismatic2DRefinement::prismatic2DRefinement
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyTopoChanger& mme
)
:
    polyMeshModifier(name, index, mme, Switch(dict.lookup("active"))),
    mesh_(mme.mesh()),
    cellsToRefine_(),
    splitPointsToUnrefine_(),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh_.nCells(), 0)
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh_.nPoints(), 0)
    ),
    refinementLevelIndicator_(0), // Must be empty before setting refinement
    level0EdgeLength_(), // Initialised in constructor body
    faceRemover_(mesh_, GREAT), // Merge boundary faces wherever possible
    maxCells_(readLabel(dict.lookup("maxCells"))),
    maxRefinementLevel_(readLabel(dict.lookup("maxRefinementLevel"))),
    edgeBasedConsistency_
    (
        dict.lookupOrDefault<Switch>("edgeBasedConsistency", true)
    ),
    nRefinementBufferLayers_
    (
        readScalar(dict.lookup("nRefinementBufferLayers"))
    ),
    nUnrefinementBufferLayers_
    (
        readScalar(dict.lookup("nUnrefinementBufferLayers"))
    )
{
    // Calculate level 0 edge length
    calcLevel0EdgeLength();

    // Check consistency between cellLevel and number of cells and pointLevel
    // and number of points in the mesh
    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorIn
        (
            "prismatic2DRefinement::prismatic2DRefinement"
            "\n("
            "\n    const word& name,"
            "\n    const dictionary& dict,"
            "\n    const label index,"
            "\n    const polyTopoChanger& mme"
            "\n)"
        )   << "Restarted from inconsistent cellLevel or pointLevel files."
            << endl
            << "Number of cells in mesh: " << mesh_.nCells()
            << " does not equal size of cellLevel: " << cellLevel_.size() << nl
            << "Number of points in mesh: " << mesh_.nPoints()
            << " does not equal size of pointLevel: " << pointLevel_.size()
            << abort(FatalError);
    }

    // Check specified number of maximum cells
    if (maxCells_ < 1)
    {
        FatalErrorIn
        (
            "prismatic2DRefinement::prismatic2DRefinement"
            "\n("
            "\n    const word& name,"
            "\n    const dictionary& dict,"
            "\n    const label index,"
            "\n    const polyTopoChanger& mme"
            "\n)"
        )   << "Specified zero or negative maxCells."
            << nl
            << "This is not allowed."
            << abort(FatalError);
    }

    // Check maximum refinement level
    if (maxRefinementLevel_ < 0)
    {
        FatalErrorIn
        (
            "prismatic2DRefinement::prismatic2DRefinement"
            "\n("
            "\n    const word& name,"
            "\n    const dictionary& dict,"
            "\n    const label index,"
            "\n    const polyTopoChanger& mme"
            "\n)"
        )   << "Negative maxRefinementLevel specified."
            << nl
            << "This is not allowed."
            << abort(FatalError);
    }

    // If the maximum refinementLevel is greater than 2 and the user insists on
    // not using point based refinement strategy, issue a warning
    if (!edgeBasedConsistency_ && maxRefinementLevel_ > 2)
    {
        WarningIn
        (
            "prismatic2DRefinement::prismatic2DRefinement"
            "\n("
            "\n    const word& name,"
            "\n    const dictionary& dict,"
            "\n    const label index,"
            "\n    const polyTopoChanger& mme"
            "\n)"
        )   << "You are not using point based consistency for dynamic"
            << " refinement."
            << nl
            << "Since you are allowing more than two maximum refinement"
            << " refinement levels, this might produce erroneous mesh due to"
            << " 8:1 point conflicts."
            << nl
            << "In order to supress this message and use point based"
            << " consistency checks, set edgeBasedConsistency to true."
            << endl;
    }

    // Check number of refinement buffer layers
    if (nRefinementBufferLayers_ < 0)
    {
        FatalErrorIn
        (
            "prismatic2DRefinement::prismatic2DRefinement"
            "\n("
            "\n    const word& name,"
            "\n    const dictionary& dict,"
            "\n    const label index,"
            "\n    const polyTopoChanger& mme"
            "\n)"
        )   << "Negative nRefinementBufferLayers specified."
            << nl
            << "This is not allowed."
            << abort(FatalError);
    }

    // Check number of unrefinement buffer layers
    if (nUnrefinementBufferLayers_ < 0)
    {
        FatalErrorIn
        (
            "prismatic2DRefinement::prismatic2DRefinement"
            "\n("
            "\n    const word& name,"
            "\n    const dictionary& dict,"
            "\n    const label index,"
            "\n    const polyTopoChanger& mme"
            "\n)"
        )   << "Negative nUnrefinementBufferLayers specified."
            << nl
            << "This is not allowed."
            << abort(FatalError);
    }

    // Check whether the number of unrefinement buffer layers is smaller than
    // number of refinement buffer layers + 2
    if (nUnrefinementBufferLayers_ < nRefinementBufferLayers_ + 2)
    {
        WarningIn
        (
            "prismatic2DRefinement::prismatic2DRefinement"
            "\n("
            "\n    const word& name,"
            "\n    const dictionary& dict,"
            "\n    const label index,"
            "\n    const polyTopoChanger& mme"
            "\n)"
        )   << "Using " << nUnrefinementBufferLayers_
            << " unrefinement buffer layers and " << nRefinementBufferLayers_
            << " refinement buffer layers."
            << nl
            << "Make sure that the number of refinement buffer layers is "
            << "at least nUnrefinementBufferLayers + 2" << nl
            << "in order to avoid problems with point level inconsistency when "
            << "refinement and unrefinement are performed in same iteration."
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::prismatic2DRefinement::~prismatic2DRefinement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::prismatic2DRefinement::setCellsToRefine
(
    const labelList& refinementCellCandidates
)
{
    if (debug)
    {
        Info<< "prismatic2DRefinement::setCellsToRefine"
            << "(const labelList& refinementCellCandidates)" << nl
            << "Setting cells to refine" << endl;
    }

    // Create a mark-up field for cells to refine
    boolList refineCell(mesh_.nCells(), false);

    // Roughly count how many cells we are going to end up with
    label roughCellCountAfterRefinement = mesh_.nCells();

    // Get cell points to count number of additional cells
    const labelListList& meshCellPoints = mesh_.cellPoints();

    // Mark initial refinement candidates for refinement only if the cell level
    // is smaller than the maximum refinement level. Note: stop marking them if
    // we exceed the rough cell count
    forAll (refinementCellCandidates, i)
    {
        // Get cell index
        const label& cellI = refinementCellCandidates[i];

        if
        (
            roughCellCountAfterRefinement < maxCells_
         && cellLevel_[cellI] < maxRefinementLevel_
        )
        {
            // Mark cell for refinement
            refineCell[cellI] = true;

            // Increment number of cells (nPoints/2 - 1 new cells per cell)
            roughCellCountAfterRefinement += meshCellPoints[cellI].size()/2 - 1;
        }
    }

    // Extend cells across faces using a specified number of refinement buffer
    // layers
    for (label i = 0; i < nRefinementBufferLayers_; ++i)
    {
        extendMarkedCellsAcrossFaces(refineCell);
    }

    // Make sure that the refinement is face consistent (2:1 consistency) and
    // point consistent (4:1 consistency) if necessary

    // Counter for additional cells to refine due to consistency in each
    // iteration and number of iterations
    label nAddCells = 0;
    label nIters = 0;
    label nTotalAddCells = 0;

    do
    {
        // Reset counter at the beginning of each iteration
        nAddCells = 0;

        if (edgeBasedConsistency_)
        {
            // Check for 4:1 edge based consistent refinement. Updates
            // cellsToRefine and returns number of cells added in this iteration
            nAddCells += edgeConsistentRefinement(refineCell);
        }

        // Check for 2:1 face based consistent refinement. Updates cellsToRefine
        // and returns number of cells added in this iteration
        nAddCells += faceConsistentRefinement(refineCell);

        // Global reduction
        reduce(nAddCells, sumOp<label>());

        // Increment number of iterations and total number of added cells
        ++nIters;
        nTotalAddCells += nAddCells;

    } while (nAddCells > 0);

    Info<< "Added " << nTotalAddCells // nTotalAddCells already reduced
        << " cells in " << returnReduce(nIters, maxOp<label>())
        << " iterations to obtain consistent refinement."
        << endl;

    // Collect all cells to refine in a dynamic list
    dynamicLabelList cellsToRefineDynamic(mesh_.nCells());

    forAll (refineCell, cellI)
    {
        if (refineCell[cellI])
        {
            // Cell marked for refinement, append it
            cellsToRefineDynamic.append(cellI);
        }
    }

    // Transfer the contents into the data member (ordinary list)
    cellsToRefine_.transfer(cellsToRefineDynamic);

    Info<< "Selected " << returnReduce(cellsToRefine_.size(), sumOp<label>())
        << " cells to refine." << endl;
}


void Foam::prismatic2DRefinement::setSplitPointsToUnrefine
(
    const labelList& unrefinementPointCandidates
)
{
    if (debug)
    {
        Info<< "prismatic2DRefinement::setSplitPointsToUnrefine"
            << "(const labelList& unrefinementPointCandidates)" << nl
            << "Setting split points to unrefine." << endl;
    }

    // Get necessary mesh data
    const label nPoints = mesh_.nPoints();
    const labelListList& meshCellPoints = mesh_.cellPoints();

    // PART 1: Mark all split points in the mesh (points that can be unrefined)
    boolList splitPointsMarkup(nPoints, false);

    // Algorithm: split point is uniquely defined as a point that:
    // 1. Has pointLevel_ > 0 (obviously),
    // 2. A point that has the same pointLevel_ as ALL of the points of its
    //    edges. In other words, for each point, we will look through all the
    //    edges of the point. For each edges, we will visit both points and
    //    check point levels. All point levels must be the same for this point
    //    candidate to be a split point. This is quite useful since there is no
    //    need to store the refinement history

    // Get necessary mesh data
    const edgeList& meshEdges = mesh_.edges();
    const labelListList& meshPointEdges = mesh_.pointEdges();

    // Loop through all points
    forAll (meshPointEdges, pointI)
    {
        // Get point level of this point
        const label& centralPointLevel = pointLevel_[pointI];

        if (centralPointLevel < 1)
        {
            // Point can't be unrefined as its level is either 0 or
            // invalid. Continue immediately
            continue;
        }

        // Flag to see whether this is a split point candidate
        bool splitPointCandidate = true;

        // Get all edge labels for this point
        const labelList& pEdges = meshPointEdges[pointI];

        // Loop through all point edges
        forAll (pEdges, i)
        {
            // Get edge index and the edge
            const label& edgeI = pEdges[i];
            const edge& curEdge = meshEdges[edgeI];

            // Loop through both points of the edge
            forAll (curEdge, j)
            {
                // Get point index
                const label& pointJ = curEdge[j];

                if (pointLevel_[pointJ] != centralPointLevel)
                {
                    // Point levels are different, this can't be a split point,
                    // set flag to false and break immediatelly
                    splitPointCandidate = false;
                    break;
                }
                // else: this is still potential split point candidate so
                //       there's nothing to do
            } // End for both points of this edge

            // Check whether this can't be a split point already and break out
            // immediately
            if (!splitPointCandidate)
            {
                break;
            }
        } // End for all point faces

        // At this point, if the flag is still true, this is a split point
        if (splitPointCandidate)
        {
            splitPointsMarkup[pointI] = true;
        }
    }

    // Note: if there is no dynamic load balancing, points at the boundary
    // cannot be split points by definition. When we implement dynamic load
    // balancing, it is possible that a split point ends up on the boundary and
    // the code should work. However, this has not been tested yet since we do
    // not have all the capability. In case there are some problems with dynamic
    // load balancing, uncomment these lines to avoid unrefining around split
    // points near processor boundaries. This might help debug the thing,
    // although I think that it should work. VV, 12/Feb/2018.
//    const label nInternalFaces = mesh_.nInternalFaces();
//    const label nFaces = mesh_.nFaces();
//
//    for (label faceI = nInternalFaces; faceI < nFaces; ++faceI)
//    {
//        // Get the face and make sure that the points are unarked
//        const face& f = meshFaces[faceI];
//
//        forAll (f, fpI)
//        {
//            splitPointsMarkup[f[fpI]] = false;
//        }
//    }

    // PART 2: Mark all unrefinement point candidates that are split points at
    // the same time (basically the intersection of split points and candidates)

    // Create markup field of split points to unrefine
    // True: this is a split point which should be unrefined
    // False: this is either not a split point or it shouldn't be unrefined
    boolList splitPointsToUnrefine(nPoints, false);

    // Loop through all unrefinement candidates
    forAll (unrefinementPointCandidates, i)
    {
        // Get point index
        const label& pointI = unrefinementPointCandidates[i];

        if (splitPointsMarkup[pointI])
        {
            // This is a split point, mark it for unrefinement
            splitPointsToUnrefine[pointI] = true;
        }
    }


    // PART 3: Make sure that we skip unrefining around split points that
    // possibly have cells around that will be refined

    // Mark cells that need to be protected (will be refined in this iteration)
    boolList protectedCell(mesh_.nCells(), false);

    // Loop through cells to refine and mark them
    forAll (cellsToRefine_, i)
    {
        protectedCell[cellsToRefine_[i]] = true;
    }

    // Extend protected cells across points using a specified number of
    // unrefinement buffer layers
    for (label i = 0; i < nUnrefinementBufferLayers_; ++i)
    {
        extendMarkedCellsAcrossPoints(protectedCell);
    }

    // Loop through all cells and if the cell should be protected, protect all
    // of its points from unrefinement
    forAll (protectedCell, cellI)
    {
        if (protectedCell[cellI])
        {
            // Get list of cell points for this protected cell
            const labelList& cPoints = meshCellPoints[cellI];

            // Loop through cell points and make sure that they are not marked
            // for unrefinement
            forAll (cPoints, j)
            {
                splitPointsToUnrefine[cPoints[j]] = false;
            }
        }
    }


    // PART 4: Ensure face consistent (2:1 constraint) and possibly edge
    // consistent (4:1 constraint) unrefinement

    // Get necessary mesh data
    const label nCells = mesh_.nCells();
    const labelListList& meshPointCells = mesh_.pointCells();

    // Count number of removed cells from unrefinement (cells that will not be
    // unrefined) in each iteration and number of iterations
    label nRemCells = 0;
    label nIters = 0;
    label nTotalRemCells = 0;

    do
    {
        // First, create cells to unrefine (all cells sharing point to unrefine)
        boolList cellsToUnrefine(nCells, false);

        // Loop through all split points to unrefine
        forAll (splitPointsToUnrefine, pointI)
        {
            if (splitPointsToUnrefine[pointI])
            {
                // This split point is marked for unrefinement, collect all of
                // its cells
                const labelList& pCells = meshPointCells[pointI];
                forAll (pCells, i)
                {
                    cellsToUnrefine[pCells[i]] = true;
                }
            }
        }

        // Reset number of removed cells from unrefinement for this iteration
        nRemCells = 0;

        if (edgeBasedConsistency_)
        {
            // Check for 4:1 edge based consistent unrefinement. Updates
            // cellsToUnrefine and returns number of removed cells from
            // unrefinement in this iteration
            nRemCells += edgeConsistentUnrefinement(cellsToUnrefine);
        }

        // Check for 2:1 face based consistent unrefinement. Updates
        // cellsToUnrefine and returns number of removed cells from unrefinement
        // in this iteration
        nRemCells += faceConsistentUnrefinement(cellsToUnrefine);

        // Global reduction
        reduce(nRemCells, sumOp<label>());

        // If we have removed at least one cell from unrefinement, we need to
        // protect its split points as well from unrefinement
        if (nRemCells > 0)
        {
            // Get point cells
            const labelListList& meshPointCells = mesh_.pointCells();

            // Loop through all split points to unrefine
            forAll (splitPointsToUnrefine, pointI)
            {
                if (splitPointsToUnrefine[pointI])
                {
                    // This is a split point for unrefinement, get the cells
                    const labelList& pCells = meshPointCells[pointI];

                    // Loop through all point cells
                    forAll (pCells, i)
                    {
                        if (!cellsToUnrefine[pCells[i]])
                        {
                            // Cell must not be refined, remove point from
                            // unrefinement as well
                            splitPointsToUnrefine[pointI] = false;
                            break;
                        }
                    }
                }
            }
        }

        // Increment number of iterations and number of total removed cells
        ++nIters;
        nTotalRemCells += nRemCells;

    } while (nRemCells > 0);

    Info<< "Removed " << nTotalRemCells // nTotalRemCells already reduced
        << " cells in " << returnReduce(nIters, maxOp<label>())
        << " iterations to obtain consistent unrefinement."
        << endl;

    // Collect all split points to unrefine in a dynamic list
    dynamicLabelList splitPointsToUnrefineDynamic(nPoints);

    forAll (splitPointsToUnrefine, pointI)
    {
        if (splitPointsToUnrefine[pointI])
        {
            // Split point marked for unrefinement, append it
            splitPointsToUnrefineDynamic.append(pointI);
        }
    }

    // Transfer the contents into the data member (ordinary list)
    splitPointsToUnrefine_.transfer(splitPointsToUnrefineDynamic);

    Info<< "Selected "
        << returnReduce(splitPointsToUnrefine_.size(), sumOp<label>())
        << " split points to unrefine." << endl;
}


bool Foam::prismatic2DRefinement::changeTopology() const
{
    if (!active())
    {
        // Modifier is inactive, skip topo change
        if (debug)
        {
            Pout<< "bool prismatic2DRefinement::changeTopology() const"
                << "for object " << name() << " : "
                << "Inactive" << endl;
        }

        return false;
    }
    else
    {
        return true;
    }
}


void Foam::prismatic2DRefinement::setRefinement(polyTopoChange& ref) const
{
    // Make sure that the point levels are updated across coupled patches before
    // setting refinement and unrefinement. Note: not sure why the sync is not
    // performed correctly if I do it in updateMesh. This is a temporary
    // solution, need to investigate in detail, but I assume something is not
    // updated yet in that case. VV, 31/Jan/2018.
    syncTools::syncPointList
    (
        mesh_,
        pointLevel_,
        maxEqOp<label>(),
        0,   // Null value
        true // Apply separation for parallel cyclics
    );

    // Set refinement and unrefinement
    setPrismatic2DRefinement(ref);
    setPrismatic2DUnrefinement(ref);

    // Clear the list of cells to refine and split points to unrefine since the
    // refinement/unrefinement instructions have been set
    cellsToRefine_.clear();
    splitPointsToUnrefine_.clear();
}


void Foam::prismatic2DRefinement::modifyMotionPoints
(
    pointField& motionPoints
) const
{
    if (debug)
    {
        Pout<< "void prismatic2DRefinement::modifyMotionPoints("
            << "pointField& motionPoints) const for object "
            << name() << " : ";
    }

    if (debug)
    {
        Pout << "No motion point adjustment" << endl;
    }
}



void Foam::prismatic2DRefinement::updateMesh(const mapPolyMesh& map)
{
    if (debug)
    {
        Info<< "prismatic2DRefinement::updateMesh(const mapPolyMesh&) "
            << " for object " << name() << " : "
            << "Updating cell and point levels."
            << endl;
    }

    // Mesh has changed topologically, we need to update cell and point levels
    // and optionally face removal object

    // Get cell map: from current mesh cells to previous mesh cells
    const labelList& cellMap = map.cellMap();

    // Create new cell level
    labelList newCellLevel(cellMap.size());

    // Loop through all new cells
    forAll (cellMap, newCellI)
    {
        // Get index of the corresponding old cell
        const label& oldCellI = cellMap[newCellI];

        if (oldCellI == -1)
        {
            // This cell is inflated (does not originate from other cell), set
            // cell level to -1
            newCellLevel[newCellI] = -1;
        }
        else
        {
            // This cell has either been added based on another cell or it
            // hasn't changed. Update new cell level according to refinement
            // level indicator and old cell level

            // Get refinement status of the old cell
            const label& refStatus = refinementLevelIndicator_[oldCellI];

            if (refStatus == UNREFINED)
            {
                // New cell has been obtained by unrefining other cells - this
                // is the remaining "master" cell. Decrement cell level
                newCellLevel[newCellI] = cellLevel_[oldCellI] - 1;
            }
            else if (refStatus == UNCHANGED)
            {
                // Cell hasn't been changed during this refinement, copy old
                // cell level
                newCellLevel[newCellI] = cellLevel_[oldCellI];
            }
            else if (refStatus == REFINED)
            {
                // Cell has been refined, increment cell level
                newCellLevel[newCellI] = cellLevel_[oldCellI] + 1;
            }
            else
            {
                FatalErrorIn
                (
                    "prismatic2DRefinement::updateMesh(const mapPolyMesh& map)"
                )   << "Invalid refinement status detected: " << refStatus << nl
                    << "Old cell index: " << oldCellI << nl
                    << "New cell index: " << newCellI << abort(FatalError);
            }
        }
    }

    // Transfer the new cell level into the data member
    cellLevel_.transfer(newCellLevel);

    // Clear out refinementLevelIndicator_ field for next refinement step
    refinementLevelIndicator_.clear();


    // Point level will be updated based on already updated cell level. Level
    // for newly added points has to be equal to the maximum cell level of
    // surrounding points. At this point, mesh is fully topologically valid so
    // it is safe to use pointCells

    // Get point map: from current mesh points to previous mesh points
    const labelList& pointMap = map.pointMap();

    // Get point cell
    const labelListList& meshPointCells = mesh_.pointCells();

    // Create new point level
    labelList newPointLevel(pointMap.size());

    // Loop through all new points
    forAll (pointMap, newPointI)
    {
        // Get index of the corresponding old point
        const label& oldPointI = pointMap[newPointI];

        if (oldPointI == -1)
        {
            // This point has been appended without any master point, use
            // surrounding cells to determine new point level

            // Get new cells surrounding this new point
            const labelList& pCells = meshPointCells[newPointI];

            // Find maximum cell level for this point
            label maxCellLevel = 0;
            forAll (pCells, i)
            {
                maxCellLevel = max(maxCellLevel, cellLevel_[pCells[i]]);
            }

            // Set new point level as the maximum of the surrounding cells
            newPointLevel[newPointI] = maxCellLevel;
        }
        else
        {
            // This point is either old point or it has been added in terms of
            // another point. New point level is equal to the old point level
            newPointLevel[newPointI] = pointLevel_[oldPointI];
        }
    }

    // Note: new point level is going to be synced at processor boundaries just
    // before the next step in setRefinement. Need to investigate why the sync
    // is not done properly if I put it here. Something is not updated yet.
    // VV, 31/Jan/2018.

    // Transfer the new point level into the data member
    pointLevel_.transfer(newPointLevel);

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Update face remover
    faceRemover_.updateMesh(map);
}


void Foam::prismatic2DRefinement::write(Ostream& os) const
{
    os  << nl << type() << nl
        << name() << nl
        << maxCells_ << nl
        << maxRefinementLevel_ << nl
        << edgeBasedConsistency_ << nl
        << nRefinementBufferLayers_ << nl
        << nUnrefinementBufferLayers_ << endl;
}


void Foam::prismatic2DRefinement::writeDict(Ostream& os) const
{
    // Write necessary data before writing dictionary
    cellLevel_.write();
    pointLevel_.write();

    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type()
        << token::END_STATEMENT << nl
        << "    maxCells " << maxCells_
        << token::END_STATEMENT << nl
        << "    maxRefinementLevel " << maxRefinementLevel_
        << token::END_STATEMENT << nl
        << "    edgeBasedConsistency " << edgeBasedConsistency_
        << token::END_STATEMENT << nl
        << "    nRefinementBufferLayers " << nRefinementBufferLayers_
        << token::END_STATEMENT << nl
        << "    nUnrefinementBufferLayers " << nUnrefinementBufferLayers_
        << token::END_STATEMENT << nl
        << "    active " << active()
        << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// ************************************************************************* //
