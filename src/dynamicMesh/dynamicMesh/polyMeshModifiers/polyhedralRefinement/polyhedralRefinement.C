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
    Generalisation of hexRef8 for polyhedral cells and refactorisation into mesh
    modifier engine.

\*---------------------------------------------------------------------------*/

#include "polyhedralRefinement.H"
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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyhedralRefinement, 0);
    addToRunTimeSelectionTable
    (
        polyMeshModifier,
        polyhedralRefinement,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::polyhedralRefinement::getAnchorLevel
(
    const label faceI
) const
{
    const face& f = mesh_.faces()[faceI];

    if (f.size() <= 3)
    {
        return pointLevel_[f[findMaxLevel(f)]];
    }
    else
    {
        const label& ownLevel = cellLevel_[mesh_.faceOwner()[faceI]];

        if (countAnchors(f, ownLevel) >= 3)
        {
            return ownLevel;
        }
        else if (countAnchors(f, ownLevel + 1) >= 3)
        {
            return ownLevel + 1;
        }
        else
        {
            return -1;
        }
    }
}


void Foam::polyhedralRefinement::calcLevel0EdgeLength()
{
    if (cellLevel_.size() != mesh_.nCells())
    {
        FatalErrorIn("scalar polyhedralRefinement::getLevel0EdgeLength() const")
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
    }

    // Get the minimum per level over all processors. Note: using minEqOp
    // because if the cells are not cubic, we end up using the smallest edge
    Pstream::listCombineGather(typEdgeLenSqr, minEqOp<scalar>());
    Pstream::listCombineScatter(typEdgeLenSqr);

    if (debug)
    {
        Pout<< "polyhedralRefinement::calcLevel0EdgeLength() :"
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
            const edge& e = meshEdges[cEdges[i]];
            const scalar edgeLenSqr = magSqr(e.vec(meshPoints));

            maxEdgeLenSqr[cLevel] = max(maxEdgeLenSqr[cLevel], edgeLenSqr);
        }
    }

    // Get maximum per level over all processors
    Pstream::listCombineGather(maxEdgeLenSqr, maxEqOp<scalar>());
    Pstream::listCombineScatter(maxEdgeLenSqr);

    if (debug)
    {
        Pout<< "polyhedralRefinement::calcLevel0EdgeLength() :"
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
        Pout<< "polyhedralRefinement::calcLevel0EdgeLength() :"
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
                Pout<< "polyhedralRefinement::calcLevel0EdgeLength() :"
                    << "Edge length: " << Foam::sqrt(lenSqr) << ", "
                    << "with edge level: " << levelI << ", "
                    << "has equivalent level0 lenght of:" << level0EdgeLength_
                    << endl;
            }

            break;
        }
    }

    if (level0EdgeLength_ == -1)
    {
        FatalErrorIn("polyhedralRefinement::calcLevel0EdgeLength()")
            << "Problem in definition of typical edge length squared: "
            << typEdgeLenSqr << abort(FatalError);
    }
}


void Foam::polyhedralRefinement::setInstance(const fileName& inst) const
{
    if (debug)
    {
        Pout<< "polyhedralRefinement::setInstance(const fileName& inst)"
            << nl
            << "Resetting file instance of refinement data to " << inst
            << endl;
    }

    cellLevel_.instance() = inst;
    pointLevel_.instance() = inst;
}


void Foam::polyhedralRefinement::extendMarkedCells(boolList& refineCell) const
{
    // Mark all faces for all marked cells
    const label nFaces = mesh_.nFaces();
    boolList refineFace(nFaces, false);

    // Get mesh cells
    const cellList& meshCells = mesh_.cells();

    // Loop through all cells
    forAll (refineCell, cellI)
    {
        if (refineCell[cellI])
        {
            // This cell is a refinement candidate, get its faces
            const cell& cFaces = meshCells[cellI];

            forAll (cFaces, i)
            {
                refineFace[cFaces[i]] = true;
            }
        }
    }

    // Snyc the face list across processor boundaries
    syncTools::syncFaceList(mesh_, refineFace, orEqOp<bool>(), false);

    // Get necessary mesh data
    const label nInternalFaces = mesh_.nInternalFaces();
    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    // Internal faces
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        if (refineFace[faceI])
        {
            // Face will be refined, mark both owner and neighbour
            refineCell[owner[faceI]] = true;
            refineCell[neighbour[faceI]] = true;
        }
    }

    // Boundary faces
    for (label faceI = nInternalFaces; faceI < nFaces; ++faceI)
    {
        if (refineFace[faceI])
        {
            // Face will be refined, mark owner
            refineCell[owner[faceI]] = true;
        }
    }
}


void Foam::polyhedralRefinement::setPolyhedralRefinement
(
    polyTopoChange& ref
) const
{
    // Note: assumes that cellsToRefine_ are set prior to the function call

    // Do nothing if there are no cells to refine
    if (cellsToRefine_.empty())
    {
        Pout<< "polyehdralRefinement::setPolyhedralRefinement(...)" << nl
            << "There are no cells selected for refinement. Returning... "
            << endl;

        return;
    }

    // Reset refinementLevelIndicator field. Note: the list is cleared in
    // updateMesh member function after updating cell and point levels
    if (refinementLevelIndicator_.empty())
    {
        // List has been resetted correctly, initialise it for this iteration
        refinementLevelIndicator_.setSize(mesh_.nCells(), UNCHANGED);
    }
    else
    {
        // List has not been resetted correctly, issue an error
        FatalErrorIn
        (
            "polyhedralRefinement::setPolyhedralRefinement(...)"
        )   << "Refinement level indicator list has not been"
            << " resetted properly." << nl
            << "Either the call to updateMesh() after performing"
            << " refinement has not been made or the call to"
            << " setPolyhedralRefinement(...) and"
            << " setPolyhedralUnrefinement(...) has not been made in"
            << " correct order." << nl
            << "Make sure to set refinement, set unrefinement and call"
            << " updateMesh after performing the topo change."
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "polyhedralRefinement::setPolyhedralRefinement(...)" << nl
            << "Allocating " << cellsToRefine_.size() << " cell midpoints."
            << endl;
    }


    // PART 1: Mark cells for refinement and add points at their cell centres

    // Get necessary mesh data
    const faceList& meshFaces = mesh_.faces();
    const cellList& meshCells = mesh_.cells();
    const vectorField& meshCellCentres = mesh_.cellCentres();

    // Mid point for refined cell (points at cell centres):
    // Not refined = -1
    // Shall be refined > -1 (label of added mid point)
    labelList cellMidPoint(mesh_.nCells(), -1);

    // Loop through cells to refine
    forAll(cellsToRefine_, i)
    {
        // Get cell idnex
        const label& cellI = cellsToRefine_[i];

        cellMidPoint[cellI] = ref.setAction
        (
            polyAddPoint
            (
                meshCellCentres[cellI], // Point to add (cell centre)
                -1,                     // Appended point: no master ID
                -1,                     // Zone for point
                true                    // Supports a cell
            )
        );
    }

    // Write out split cells as a cell set for debug
    if (debug)
    {
        // Note: cellSet is actually a hash table of labels
        cellSet splitCells(mesh_, "splitCells", cellsToRefine_.size());

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] > -1)
            {
                // Cell is marked for refinement, insert into cellSet
                splitCells.insert(cellI);
            }
        }

        Pout<< "polyhedralRefinement::setPolyhedralRefinement(...)" << nl
            << "Writing " << splitCells.size()
            << " cells to split to cellSet " << splitCells.objectPath()
            << endl;

        splitCells.write();
    }


    // PART 2: Mark edges for refinement and add points to edge centres

    if (debug)
    {
        Pout<< "polyhedralRefinement::setPolyhedralRefinement(...)" << nl
            << "Allocating edge midpoints."
            << endl;
    }

    // Refined edges are defined by having both their point levels lower than
    // the cell level, i.e. if any cell that gets split uses this edge, the edge
    // needs to be split as well

    // Get necessary mesh data
    const labelListList& meshCellEdges = mesh_.cellEdges();
    const edgeList& meshEdges = mesh_.edges();

    // Mid points for refined edge:
    // No need to split edge = -1
    // Label of introduced mid point > -1
    labelList edgeMidPoint(mesh_.nEdges(), -1);

    // Note: Loop over all cells instead of all edges
    forAll(cellMidPoint, cellI)
    {
        // Point is marked for refinement, proceed to look at the edges
        if (cellMidPoint[cellI] > -1)
        {
            // Get edges of this cell
            const labelList& cEdges = meshCellEdges[cellI];

            forAll(cEdges, i)
            {
                // Get edge index and edge
                const label& edgeI = cEdges[i];
                const edge& e = meshEdges[edgeI];

                if
                (
                    pointLevel_[e[0]] <= cellLevel_[cellI]
                 && pointLevel_[e[1]] <= cellLevel_[cellI]
                )
                {
                    // Point levels of both edge points are <= cell level, mark
                    // the edge for splitting
                    edgeMidPoint[edgeI] = 12345;
                }
            }
        }
    }

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

    // Introduce edge points

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

        Pout<< "polyhedralRefinement::setPolyhedralRefinement(...)" << nl
            << "Writing centres of edges to split to file " << str.name()
            << endl;
    }


    // PART 3: Calculate face level (after selected cells splitting)

    if (debug)
    {
        Pout<< "polyhedralRefinement::setPolyhedralRefinement" << nl
            << "Allocating face midpoints."
            << endl;
    }

    // Face anchor level. There are guaranteed at least 3 points with level
    // <= anchorLevel. These are the corner points.
    labelList faceAnchorLevel(mesh_.nFaces());

    for (label faceI = 0; faceI < mesh_.nFaces(); ++faceI)
    {
        faceAnchorLevel[faceI] = getAnchorLevel(faceI);
    }

    // Mid points for refined face (points at face centres):
    // Not refined = -1
    // Shall be refined > -1 (label of added mid point)
    labelList faceMidPoint(mesh_.nFaces(), -1);

    // Get necessary mesh data
    const labelList& meshFaceOwner = mesh_.faceOwner();
    const labelList& meshFaceNeighbour = mesh_.faceNeighbour();

    const label nFaces = mesh_.nFaces();
    const label nInternalFaces = mesh_.nInternalFaces();

    // Internal faces: look at cells on both sides. Uniquely determined since
    // the face itself is guaranteed to be same level as most refined neighbour
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        // Note: no need to check whether the face has valid anchor level since
        // all faces can be split
        const label& own = meshFaceOwner[faceI];
        const label& ownLevel = cellLevel_[own];
        const label newOwnLevel = ownLevel + (cellMidPoint[own] > -1 ? 1 : 0);

        const label& nei = meshFaceNeighbour[faceI];
        const label& neiLevel = cellLevel_[nei];
        const label newNeiLevel = neiLevel + (cellMidPoint[nei] > -1 ? 1 : 0);

        if
        (
            newOwnLevel > faceAnchorLevel[faceI]
         || newNeiLevel > faceAnchorLevel[faceI]
        )
        {
            // New level is higher than the face anchor level, mark for
            // splitting
            faceMidPoint[faceI] = 12345;
        }
    }

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
            const label newOwnLevel =
                ownLevel + (cellMidPoint[own] > -1 ? 1 : 0);

            newNeiLevel[i] = newOwnLevel;
        }

        // Swap the list which now contains data from the other side
        syncTools::swapBoundaryFaceList(mesh_, newNeiLevel, false);

        forAll(newNeiLevel, i)
        {
            // Get face index
            const label faceI = i + nInternalFaces;

            // Note: no need to check whether the face has valid anchor level
            // since all faces can be split
            const label& own = meshFaceOwner[faceI];
            const label& ownLevel = cellLevel_[own];
            const label newOwnLevel =
                ownLevel + (cellMidPoint[own] > -1 ? 1 : 0);

            if
            (
                newOwnLevel > faceAnchorLevel[faceI]
             || newNeiLevel[i] > faceAnchorLevel[faceI]
            )
            {
                // New level is higher than the face anchor level, mark for
                // splitting
                faceMidPoint[faceI] = 12345;
            }
        }
    } // End memory management for syncing owner/neighbour face levels

    // Note: synronisation of faceMidPoints across coupled patches is not
    // necessary since we have exchanged the neighbour data above using
    // swapBoundaryFaceList, thus the faceMidPoint has to be the same on both
    // sides. VV, 5/Jan/2018.
//    // Synchronize faceMidPoint across coupled patches
//    syncTools::syncFaceList
//    (
//        mesh_,
//        faceMidPoint,
//        maxEqOp<label>(),
//        false
//    );


    // Introduce face points

    // Get face centres
    const vectorField& meshFaceCentres = mesh_.faceCentres();

    // Memory management
    {
        // Phase 1: determine mid points and sync. Note: the same procedure has
        // been used for syncing edge mid points

        // Allocate storage for boundary face points
        pointField bFaceMids
        (
            nFaces - nInternalFaces,
            point(-GREAT, -GREAT, -GREAT)
        );

        // Loop through boundary face mids
        forAll(bFaceMids, i)
        {
            // Get face index
            const label faceI = i + nInternalFaces;

            if (faceMidPoint[faceI] > -1)
            {
                // This is a valid face mid, get the face centre
                bFaceMids[i] = meshFaceCentres[faceI];
            }
        }

        // Sync across coupled boundaries. Note: uses maximum of the components
        // of the vector. This is completely arbitrary but it doesn't matter as
        // long as we have same points on each side. VV, 5/Jan/2018.
        syncTools::syncBoundaryFaceList
        (
            mesh_,
            bFaceMids,
            maxEqOp<vector>(),
            true               // apply separation
        );

        // Loop through faces
        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] > -1)
            {
                // Face marked to be split. Add the point at face centre and
                // replace faceMidPoint with actual point label

                faceMidPoint[faceI] = ref.setAction
                (
                    polyAddPoint
                    (
                        (
                            faceI < nInternalFaces
                          ? meshFaceCentres[faceI]
                          : bFaceMids[faceI - nInternalFaces]
                        ),    // Point
                        -1,   // Appended point, no master ID
                        -1,   // Zone for point
                        true  // Supports a cell
                    )
                );
            }
        }
    } // End memory management for syncing boundary data and adding face mids

    // Write out split faces as a face set for debugging
    if (debug)
    {
        // Create a faceSet with 3*cell sizes to prevent excessive resizing
        faceSet splitFaces(mesh_, "splitFaces", 3*cellsToRefine_.size());

        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] > -1)
            {
                splitFaces.insert(faceI);
            }
        }

        Pout<< "polyhedralRefinement::setPolyhedralRefinement(...)" << nl
            << "Writing " << splitFaces.size()
            << " faces to split to faceSet " << splitFaces.objectPath()
            << endl;

        splitFaces.write();
    }


    // Now we have all the information we need to perform the refinement and we
    // no longer need to refer to cellsToRefine_. The information is:
    // - cellMidPoint >= 0 : cell needs to be split
    // - faceMidPoint >= 0 : face needs to be split
    // - edgeMidPoint >= 0 : edge needs to be split


    // PART 4: Get corner and anchor points for all cells

    if (debug)
    {
        Pout<< "polyhedralRefinement::setPolyhedralRefinement(...)" << nl
            << "Finding cell anchorPoints" << endl;
    }

    // Get anchor points for each cell: points that have the same or lower
    // refinement level as the cell
    List<dynamicLabelList> cellAnchorPointsDynamic(mesh_.nCells());

    // Loop through all cells
    forAll(cellMidPoint, cellI)
    {
        if (cellMidPoint[cellI] > -1)
        {
            // This cell shall be cut. Set capacity to 16 to prevent excessive
            // resizing
            cellAnchorPointsDynamic[cellI].setCapacity(16);
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
                cellMidPoint[cellI] > -1
             && pointLevel_[pointI] <= cellLevel_[cellI]
            )
            {
                // This point cells is marked for refinement and its point level
                // is smaller or equal to cell level, append the point
                cellAnchorPointsDynamic[cellI].append(pointI);
            }
        }
    }

    // Loop through all cells and check whether at least 4 anchor points
    // have been found (minimum requirement for a tet cell)

    // Get cell points for error output
    const labelListList& meshCellPoints = mesh_.cellPoints();

    forAll(cellMidPoint, cellI)
    {
        if (cellMidPoint[cellI] > -1)
        {
            // Cell selected for refinement
            if (cellAnchorPointsDynamic[cellI].size() < 4)
            {
                // Cell has less than 4 anchor points. Issue an error and report
                // cell points
                const labelList& cPoints = meshCellPoints[cellI];

                FatalErrorIn
                (
                    "polyhedralRefinement::setPolyhedralRefinement(...)"
                )   << "Cell " << cellI
                    << " of level " << cellLevel_[cellI]
                    << " does not seem to have enough points of "
                    << " lower level" << endl
                    << "cellPoints:" << cPoints << endl
                    << "pointLevels:"
                    << IndirectList<label>(pointLevel_, cPoints)() << endl
                    << abort(FatalError);
            }
        }
    }

    // Collect cellAnchorPoints into a List<labelList> instead of
    // List<dynamicList>
    labelListList cellAnchorPoints(mesh_.nCells());

    forAll(cellAnchorPointsDynamic, cellI)
    {
        // Tranfer the dynamic list for each cell into an ordinary list
        cellAnchorPoints[cellI].transfer(cellAnchorPointsDynamic[cellI]);
    }

    // PART 5: Add the cells

    if (debug)
    {
        Pout<< "polyhedralRefinement::setPolyhedralRefinement(...)" << nl
            << " Adding cells."
            << endl;
    }

    // We should have exactly n new cells per each split cell, where n is the
    // number of anchor points in a cell
    labelListList cellAddedCells(mesh_.nCells());

    // Get cell zone mesh
    const cellZoneMesh& cellZones = mesh_.cellZones();

    forAll(cellAnchorPoints, cellI)
    {
        // Check whether this is a split cell
        if (cellMidPoint[cellI] > -1)
        {
            // Get cell anchors
            const labelList& cAnchors = cellAnchorPoints[cellI];

            // Set the total number of added cells to number of anchors
            labelList& cAdded = cellAddedCells[cellI];
            cAdded.setSize(cAnchors.size());

            // Original cell has index 0
            cAdded[0] = cellI;

            // Update refinement level indicator field to 1 since this original
            // cell will be refined
            refinementLevelIndicator_[cellI] = REFINED;

            // Add other cells
            for (label i = 1; i < cAdded.size(); ++i)
            {
                cAdded[i] = ref.setAction
                (
                    polyAddCell
                    (
                        -1,                         // Master point
                        -1,                         // Master edge
                        -1,                         // Master face
                        cellI,                      // Master cell
                        cellZones.whichZone(cellI)  // Zone for cell
                    )
                );
            }
        }
    }


    // PART 6: Adding faces

    // 6.1. Existing faces that get split (into n faces where n is the number of
    //      points or edges)
    // 6.2. Existing faces that do not get split but only edges get split
    // 6.3. Existing faces that do not get split but get new owner/neighbour
    // 6.4. New internal faces inside split cells.

    if (debug)
    {
        Pout<< "polyhedralRefinement::setPolyhedralRefinement(...)" << nl
            << " Marking faces to be handled"
            << endl;
    }

    // Get all faces to split:
    // a) All faces of a cell being split
    // b) All faces that are being split
    // c) Both faces of an edge that is being split
    boolList facesToSplit(mesh_.nFaces(), false);

    // Get edge faces
    const labelListList& meshEdgeFaces = mesh_.edgeFaces();

    // a) All faces of a cell that is being split
    forAll(cellMidPoint, cellI)
    {
        if (cellMidPoint[cellI] > -1)
        {
            const cell& cFaces = meshCells[cellI];

            forAll(cFaces, i)
            {
                facesToSplit[cFaces[i]] = true;
            }
        }
    }

    // b) All faces that are being split
    forAll(faceMidPoint, faceI)
    {
        if (faceMidPoint[faceI] > -1)
        {
            facesToSplit[faceI] = true;
        }
    }

    // c) Both faces of an edge that are being split
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


    // PART 6.1. Add/modify faces for each face being split

    if (debug)
    {
        Pout<< "polyhedralRefinement::setPolyhedralRefinement(...)" << endl
            << "Splitting faces..." << endl;
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
            forAll(f, fp)
            {
                const label& pointI = f[fp];

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
                        fp,
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
                        fp,
                        faceVerts
                    );

                    // Transfer dynamic list to a face (ordinary list)
                    newFace.transfer(faceVerts);
                    faceVerts.clear();

                    // Set new owner/neighbour indices based on split cells
                    label own, nei;
                    setNewFaceNeighbours
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        faceI,
                        pointI,           // Anchor point

                        own,
                        nei
                    );


                    if (debug)
                    {
                        // Get mesh cell centres
                        const vectorField& meshCellCentres =
                            mesh_.cellCentres();

                        if (mesh_.isInternalFace(faceI))
                        {
                            const label oldOwn = meshFaceOwner[faceI];
                            const label oldNei = meshFaceNeighbour[faceI];

                            // Print info only with deep debug level
                            if (debug > 1)
                            {
                                Pout<< "Split internal face: " << faceI
                                    << ", verts: " << f
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
                            const label oldOwn = meshFaceOwner[faceI];

                            // Print info only with deep debug level
                            if (debug > 1)
                            {
                                Pout<< "Split boundary face: " << faceI
                                    << ", verts: " << f
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
                    } // End debug


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
        }
    }


    // PART 6.2. Modify faces that do not get split but have edges that are
    // being split

    if (debug)
    {
        Pout<< "polyhedralRefinement::setPolyhedralRefinement(...)" << nl
            << "Modifying faces with split edges..."
            << endl;
    }

    // Get face edges
    const labelListList& meshFaceEdges = mesh_.faceEdges();

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

                // Check whether this is not a face that has been marked for
                // splitting and that the face has not been handled yet. The
                // second check is necessary since we go through edge faces
                // instead of just faces
                if (faceMidPoint[faceI] < 0 && facesToSplit[faceI])
                {
                    // This is unsplit face that has not been handled

                    // Get face and face edges
                    const face& f = meshFaces[faceI];
                    const labelList& fEdges = meshFaceEdges[faceI];

                    // Create a dynamic list containing new face vertices
                    dynamicLabelList newFaceVerts(f.size());

                    // Append all original points and all edge mid points
                    forAll(f, fp)
                    {
                        newFaceVerts.append(f[fp]);

                        const label edgeI = fEdges[fp];

                        if (edgeMidPoint[edgeI] > -1)
                        {
                            newFaceVerts.append(edgeMidPoint[edgeI]);
                        }
                    }

                    face newFace;
                    newFace.transfer(newFaceVerts.shrink());


                    // The point with the lowest level should be an anchor
                    // point of the neighbouring cells.
                    const label anchorFp = findMinLevel(f);

                    label own, nei;
                    setNewFaceNeighbours
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        faceI,
                        f[anchorFp],      // Anchor point

                        own,
                        nei
                    );


                    if (debug)
                    {
                        // Get mesh cell centres
                        const vectorField& meshCellCentres =
                            mesh_.cellCentres();

                        if (mesh_.isInternalFace(faceI))
                        {
                            const label oldOwn = meshFaceOwner[faceI];
                            const label oldNei = meshFaceNeighbour[faceI];

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
                            const label oldOwn = meshFaceOwner[faceI];

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
                    } // End debug

                    // Modify the face
                    modifyFace(ref, faceI, newFace, own, nei);

                    // Mark face as handled
                    facesToSplit[faceI] = false;

                }// End if unsplit, unhandled face
            } // End for all edge faces
        } // End if edge has been cut
    } // End for all edges


    // PART 6.3.: Modify faces that do not get split but whose owner/neighbour
    // change due to splitting

    if (debug)
    {
        Pout<< "polyhedralRefinement::setPolyhedralRefinement(...)" << nl
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
            label anchorFp = findMaxLevel(f);

            label own, nei;
            setNewFaceNeighbours
            (
                cellAnchorPoints,
                cellAddedCells,
                faceI,
                f[anchorFp],      // Anchor point

                own,
                nei
            );

            // Modify the face, changing owner and neighbour
            modifyFace(ref, faceI, f, own, nei);

            // Mark face as handled
            facesToSplit[faceI] = false;
        }
    }


    // PART 6.4. Addd new internal faces inside split cells

    // We have to find the splitting points between the anchor points. But the
    // edges between the anchor points might have been split (into two, three or
    // four edges)

    if (debug)
    {
        Pout<< "polyhedralRefinement::setPolyhedralRefinement(...)" << nl
            << " Adding new internal faces for split cells..."
            << endl;
    }

    // Loop through cells
    forAll(cellMidPoint, cellI)
    {
        if (cellMidPoint[cellI] > -1)
        {
            // Cell has been split, create internal faces
            createInternalFaces
            (
                cellAnchorPoints,
                cellAddedCells,
                cellMidPoint,
                faceMidPoint,
                faceAnchorLevel,
                edgeMidPoint,
                cellI,
                ref
            );
        }
    }

    // Debug: check minimum point index of added points, needs to be equal to
    // number of points in the original mesh
    if (debug)
    {
        label minPointI = labelMax;
        label maxPointI = labelMin;

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] > -1)
            {
                minPointI = min(minPointI, cellMidPoint[cellI]);
                maxPointI = max(maxPointI, cellMidPoint[cellI]);
            }
        }
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
            FatalErrorIn("polyhedralRefinement::setPolyhedralRefinement(...)")
                << "Added point labels not consecutive to existing mesh points."
                << nl
                << "mesh_.nPoints():" << mesh_.nPoints()
                << " minPointI: " << minPointI
                << " maxPointI: " << maxPointI
                << abort(FatalError);
        }
    }
}


void Foam::polyhedralRefinement::setPolyhedralUnrefinement
(
    polyTopoChange& ref
) const
{
    // Resize refinementLevelIndicator field if necessary
    if (refinementLevelIndicator_.empty())
    {
        // The list is empty: meaning that we do not have any cells to refine in
        // this iteration. Resize the list and mark all cells as unchanged
        refinementLevelIndicator_.setSize(mesh_.nCells(), UNCHANGED);
    }
    // else we have some cells to refine and the list has already been set in
    // setPolyhedralRefinement member function

    // Get point cells necessary for debug and face removal
    const labelListList& meshPointCells = mesh_.pointCells();

    if (debug)
    {
        Pout<< "polyhedralRefinement::setPolyhedralUnrefinement"
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
                    "polyhedralRefinement::setPolyhedralUnrefinement"
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

        Pout<< "polyhedralRefinement::setPolyhedralUnrefinement"
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
        // Collect split faces in the hast set, guess size to prevent excessive
        // resizing
        labelHashSet splitFaces(12*splitPointsToUnrefine_.size());

        // Get point faces
        const labelListList& meshPointFaces = mesh_.pointFaces();

        forAll(splitPointsToUnrefine_, i)
        {
            // Loop through all faces of this point insert face index
            const labelList& pFaces = meshPointFaces[splitPointsToUnrefine_[i]];

            forAll(pFaces, j)
            {
                splitFaces.insert(pFaces[j]);
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
                "polyhedralRefinement::setPolyhedralUnrefinement"
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


Foam::label Foam::polyhedralRefinement::addFace
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


Foam::label Foam::polyhedralRefinement::addInternalFace
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


void Foam::polyhedralRefinement::modifyFace
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


void Foam::polyhedralRefinement::createInternalFaces
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
    const cell& meshCells = mesh_.cells()[cellI];
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
    forAll(meshCells, i)
    {
        // Get face index
        const label& faceI = meshCells[i];

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
                            "polyhedralRefinement::createInternalFaces(...)"
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
                            "polyhedralRefinement::createInternalFaces(...)"
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


Foam::label Foam::polyhedralRefinement::getAnchorCell
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

        FatalErrorIn("polyhedralRefinement::getAnchorCell(...)")
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


void Foam::polyhedralRefinement::setNewFaceNeighbours
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label faceI,
    const label pointI,

    label& own,
    label& nei
) const
{
    // Get anchor cell for this anchor point on owner side
    own = getAnchorCell
    (
        cellAnchorPoints,
        cellAddedCells,
        mesh_.faceOwner()[faceI],
        faceI,
        pointI
    );

    if (mesh_.isInternalFace(faceI))
    {
        // Get anchor cell for this anchor point on neighbour side
        nei = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            mesh_.faceNeighbour()[faceI],
            faceI,
            pointI
        );
    }
    else
    {
        // Boundary face: set neighbour to -1
        nei = -1;
    }
}


void Foam::polyhedralRefinement::walkFaceToMid
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


void Foam::polyhedralRefinement::walkFaceFromMid
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


Foam::label Foam::polyhedralRefinement::findMinLevel(const labelList& f) const
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


Foam::label Foam::polyhedralRefinement::findMaxLevel(const labelList& f) const
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


Foam::label Foam::polyhedralRefinement::countAnchors
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


Foam::label Foam::polyhedralRefinement::findLevel
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
            FatalErrorIn("polyhedralRefinement::findLevel(...)")
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

    FatalErrorIn("polyhedralRefinement::findLevel(...)")
        << "face:" << f
        << " level:" << IndirectList<label>(pointLevel_, f)()
        << " startFp:" << startFp
        << " wantedLevel:" << wantedLevel
        << abort(FatalError);

    return -1;
}


Foam::label Foam::polyhedralRefinement::storeMidPointInfo
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


void Foam::polyhedralRefinement::insertEdgeSplit
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


void Foam::polyhedralRefinement::checkInternalOrientation
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
            "polyhedralRefinement::checkInternalOrientation(...)"
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
            "polyhedralRefinement::checkInternalOrientation(...)"
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


void Foam::polyhedralRefinement::checkBoundaryOrientation
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
            "polyhedralRefinement::checkBoundaryOrientation(...)"
        )
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
            "polyhedralRefinement::checkBoundaryOrientation(...)"
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


Foam::label Foam::polyhedralRefinement::faceConsistentRefinement
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

        if (neiLevel[i] > (curOwnLevel + 1))
        {
            // Neighbour level is higher than owner level + 1, owner must be
            // marked for refinement
            cellsToRefine[own] = true;
            ++nAddCells;
        }
        // Note: other possibility (that owner level is higher than neighbour
        // level + 1) is taken into account on the other side automatically
    }

    // Return number of added cells
    return nAddCells;
}


Foam::label Foam::polyhedralRefinement::pointConsistentRefinement
(
    boolList& cellsToRefine
) const
{
    // Count number of cells that will be added
    label nAddCells = 0;

    // Collect all points from cells to refine. Assume that 10% of mesh poitns
    // are going to be affected to prevent excessive resizing
    labelHashSet pointsToConsider(mesh_.nPoints()/10);

    // Get cell points
    const labelListList& meshCellPoints = mesh_.cellPoints();

    // Collect all points
    forAll (meshCellPoints, cellI)
    {
        if (cellsToRefine[cellI])
        {
            // Get current cell points and insert them into hash set
            const labelList& curPoints = meshCellPoints[cellI];
            forAll (curPoints, pointI)
            {
                pointsToConsider.insert(curPoints[pointI]);
            }
        }
    }

    // Maximum surrounding cell refinement level for each point
    labelList maxRefLevel(mesh_.nPoints(), 0);

    // Get point cells
    const labelListList& meshPointCells = mesh_.pointCells();

    // Loop through all points and collect maximum point level for each point
    forAllConstIter (labelHashSet, pointsToConsider, iter)
    {
        // Get point index
        const label& pointI = iter.key();

        // Get the cells for this point
        const labelList& curCells = meshPointCells[pointI];

        // Find maximum refinement level for this points
        forAll (curCells, cellI)
        {
            // Get cell index and cell level
            const label& curCellI = curCells[cellI];
            const label curLevel =
                cellsToRefine[curCellI]
              ? cellLevel_[curCellI] + 1
              : cellLevel_[curCellI];

            // Set the maximum point refinement level
            if (curLevel > maxRefLevel[pointI])
            {
                maxRefLevel[pointI] = curLevel;
            }
        }
    }

    // Sync maximum refinement level across coupled boundaries
    syncTools::syncPointList
    (
        mesh_,
        maxRefLevel,
        maxEqOp<label>(),
        0,   // Null value
        true // Apply separation for parallel cyclics
    );

    // Now that the levels are synced, go through considered points and add
    // cells to refine
    forAllConstIter (labelHashSet, pointsToConsider, iter)
    {
        // Get point index
        const label& pointI = iter.key();

        // Get the cells for this point
        const labelList& curCells = meshPointCells[pointI];

        // Loop through these point cells and set cells for refinement which
        // would end up having refinement level smaller than maximum level - 1
        forAll (curCells, cellI)
        {
            // Get cell index, reference to refinement flag and cell level
            const label& curCellI = curCells[cellI];
            bool& willBeRefined = cellsToRefine[curCellI];
            const label curLevel =
                willBeRefined
              ? cellLevel_[curCellI] + 1
              : cellLevel_[curCellI];

            if (curLevel < maxRefLevel[pointI] - 1)
            {
                if (!willBeRefined)
                {
                    // Cell has not been marked for refinement, mark the cell for
                    // refinement and increment the counter
                    willBeRefined = true;
                    ++nAddCells;
                }
                else
                {
                    FatalErrorIn
                    (
                        "label polyhedralRefinement::pointConsistentRefinement"
                        "(boolList cellsToRefine) const"
                    )   << "Cell is marked for refinement, but the 4:1 point"
                        << " consistency cannot be ensured." << nl
                        << "Something went wrong before this step."
                        << endl;
                }
            }
        }
    }

    return nAddCells;
}


Foam::label Foam::polyhedralRefinement::faceConsistentUnrefinement
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

        if (ownLevel < (neiLevel -1))
        {
            // Owner level is smaller than neighbour level - 1, we must not
            // unrefine owner

            // Check whether the cell has not been marked for unrefinement
            if (!cellsToUnrefine[own])
            {
                FatalErrorIn
                (
                    "label polyhedralRefinement::faceConsistentUnrefinement"
                    "(boolList& cellsToUnrefine)"
                )   << "Cell not marked for unrefinement, indicating a"
                    << " previous unnoticed problem with unrefinement."
                    << nl
                    << "Owner: " << own << ", neighbour: " << nei
                    << nl
                    << "Owner level: " << ownLevel
                    << ", neighbour level: " << neiLevel
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
                    "label polyhedralRefinement::faceConsistentUnrefinement"
                    "(boolList& cellsToUnrefine)"
                )   << "Cell not marked for unrefinement, indicating a"
                    << " previous unnoticed problem with unrefinement."
                    << nl
                    << "Owner: " << own << ", neighbour: " << nei
                    << nl
                    << "Owner level: " << ownLevel
                    << ", neighbour level: " << neiLevel
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

        if (curOwnLevel < (neiLevel[i] - 1))
        {
            // Owner level is smaller than neighbour level - 1, we must not
            // unrefine owner

            // Check whether the cell has not been marked for unrefinement
            if (!cellsToUnrefine[own])
            {
                FatalErrorIn
                (
                    "label polyhedralRefinement::faceConsistentUnrefinement"
                    "(boolList& cellsToUnrefine)"
                )   << "Boundary cell not marked for unrefinement, indicating a"
                    << " previous unnoticed problem with unrefinement."
                    << nl
                    << "Owner: " << own
                    << nl
                    << "Owner level: " << curOwnLevel
                    << ", neighbour level: " << neiLevel[i]
                    << abort(FatalError);
            }

            cellsToUnrefine[own] = false;
            ++nRemCells;
        }

        // Note: other possibility (that neighbour level is smaller than owner
        // level - 1) is taken into account on the other side automatically
    }

    // Return number of local cells removed from unrefinement
    return nRemCells;
}


Foam::label Foam::polyhedralRefinement::pointConsistentUnrefinement
(
    const boolList& splitPointsToUnrefine,
    boolList& cellsToUnrefine
) const
{
    // Count number of cells removed from unrefinement
    label nRemCells = 0;

    // Get necessary mesh data
    const label nPoints = mesh_.nPoints();
    const labelListList& meshPointCells = mesh_.pointCells();

    // Minimum cell refinement level for each point. Note: initialise with
    // labelMax
    labelList minRefLevel(nPoints, labelMax);

    // Loop through all points and collect minimum point level for each point
    forAll (splitPointsToUnrefine, pointI)
    {
        // Check whether the point is marked for unrefinement
        if (splitPointsToUnrefine[pointI])
        {
            // Get cells for this point
            const labelList& curCells = meshPointCells[pointI];

            // Find minimum refinement level for this points
            forAll (curCells, cellI)
            {
                // Get current cell and its new level
                const label curCellI = curCells[cellI];
                const label curLevel =
                    cellsToUnrefine[curCellI]
                  ? cellLevel_[curCellI] - 1
                  : cellLevel_[curCellI];

                // If the current level is smaller than current minimum
                // refinement level, update the minimum refinement level
                if (curLevel < minRefLevel[pointI])
                {
                    minRefLevel[pointI] = curLevel;
                }
            }
        }
    }

    // Sync minimum refinement level across coupled boundaries
    syncTools::syncPointList
    (
        mesh_,
        minRefLevel,
        minEqOp<label>(),
        0,   // Null value
        true // Apply separation for parallel cyclics
    );

    // Now that the levels are synced, go through split point candidates and add
    // cells to unrefine
    forAll (splitPointsToUnrefine, pointI)
    {
        // Check whether the point is marked for unrefinement
        if (splitPointsToUnrefine[pointI])
        {
            // Get the cells for this point
            const labelList& curCells = meshPointCells[pointI];

            // Loop through these point cells and set cells for unrefinement which
            // would end up having refinement level greater than level + 1
            forAll (curCells, cellI)
            {
                // Get current cell index and level
                const label curCellI = curCells[cellI];
                const label curLevel =
                    cellsToUnrefine[curCellI]
                  ? cellLevel_[curCellI] - 1
                  : cellLevel_[curCellI];

                if (curLevel > minRefLevel[pointI] + 1)
                {
                    // Check whether the cell has not been marked for unrefinement
                    if (!cellsToUnrefine[curCellI])
                    {
                        // Set the cell for unrefinement and increment the
                        // counter
                        cellsToUnrefine[curCellI] = true;
                        ++nRemCells;
                    }
                    else
                    {
                        FatalErrorIn
                        (
                            "label"
                            "polyhedralRefinement::pointConsistentUnrefinement"
                            "\n("
                            "\n    const boolList& splitPointsToUnrefine,"
                            "\n    boolList& cellsToUnrefine"
                            "\n) const"
                        )   << "Cell is marked for unrefinement, but the 4:1 point"
                            << " consistency cannot be ensured." << nl
                            << "Something went wrong before this step."
                            << endl;
                    }
                }
            }
        }
    }

    return nRemCells;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyhedralRefinement::polyhedralRefinement
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
    pointBasedRefinement_
    (
        dict.lookupOrDefault<Switch>("pointBasedRefinement", true)
    ),
    nBufferLayers_(readScalar(dict.lookup("nBufferLayers")))
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
            "polyhedralRefinement::polyhedralRefinement"
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
            "polyhedralRefinement::polyhedralRefinement"
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
            "polyhedralRefinement::polyhedralRefinement"
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
    if (!pointBasedRefinement_ && maxRefinementLevel_ > 2)
    {
        WarningIn
        (
            "polyhedralRefinement::polyhedralRefinement"
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
            << " consistency checks, set pointBasedRefinement to true."
            << endl;
    }

    // Check number of buffer layers
    if (nBufferLayers_ < 0)
    {
        FatalErrorIn
        (
            "polyhedralRefinement::polyhedralRefinement"
            "\n("
            "\n    const word& name,"
            "\n    const dictionary& dict,"
            "\n    const label index,"
            "\n    const polyTopoChanger& mme"
            "\n)"
        )   << "Negative nBufferLayers specified."
            << nl
            << "This is not allowed."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyhedralRefinement::~polyhedralRefinement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyhedralRefinement::setCellsToRefine
(
    const labelList& refinementCellCandidates
)
{
    if (debug)
    {
        Info<< "polyhedralRefinement::setCellsToRefine"
            << "(const labelList& refinementCellCandidates)" << nl
            << "Setting cells to refine" << endl;
    }

    if (refinementCellCandidates.empty())
    {
        // Set cells to refine to empty list and return
        cellsToRefine_.clear();
        return;
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

            // Increment number of cells (nPoints - 1 new cells per cell)
            roughCellCountAfterRefinement += meshCellPoints[cellI].size() - 1;
        }
    }

    // Extend cells using a specified number of buffer layers
    for (label i = 0; i < nBufferLayers_; ++i)
    {
        extendMarkedCells(refineCell);
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

        if (pointBasedRefinement_)
        {
            // Check for 4:1 point based consistent refinement. Updates
            // cellsToRefine and returns number of cells added in this iteration
            nAddCells += pointConsistentRefinement(refineCell);
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

    Info<< "polyhedralRefinement::setCellsToRefine"
        << "(const labelList& refinementCellCandidates)" << nl
        << "Selected " << returnReduce(cellsToRefine_.size(), sumOp<label>())
        << " cells to refine." << endl;
}


void Foam::polyhedralRefinement::setSplitPointsToUnrefine
(
    const labelList& unrefinementPointCandidates
)
{
    if (debug)
    {
        Info<< "polyhedralRefinement::setSplitPointsToUnrefine"
            << "(const labelList& unrefinementPointCandidates)" << nl
            << "Setting split points to unrefine." << endl;
    }

    if (unrefinementPointCandidates.empty())
    {
        // Set split points to unrefine to empty list and return
        splitPointsToUnrefine_.clear();
        return;
    }

    // Get necessary mesh data
    const label nPoints = mesh_.nPoints();
    const labelListList& meshCellPoints = mesh_.cellPoints();

    // PART 1: Mark all split points in the mesh (points that can be unrefined)
    boolList splitPointsMarkup(nPoints, false);

    // Algorithm: split point is uniquely defined as a point that:
    // 1. Has pointLevel_ > 0 (obviously),
    // 2. A point that has the same pointLevel_ as ALL of the points of its
    //    faces. In other words, for each point, we will look through all the
    //    faes of the point. For each of the face, we will visit points and
    //    check the point level of all of these points. All point levels must be
    //    the same for this point candidate to be split point. This is quite
    //    useful since there is no need to store the refinement history

    // Get necessary mesh data
    const faceList& meshFaces = mesh_.faces();
    const labelListList& meshPointFaces = mesh_.pointFaces();

    // Loop through all points
    forAll (meshPointFaces, pointI)
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

        // Get face labels for this point
        const labelList& pFaces = meshPointFaces[pointI];

        // Loop through all point faces
        forAll (pFaces, i)
        {
            // Get face index and the face
            const label& faceI = pFaces[i];
            const face& curFace = meshFaces[faceI];

            // Loop through points of the face
            forAll (curFace, j)
            {
                // Get point index
                const label& pointJ = curFace[j];

                if (pointLevel_[pointJ] != centralPointLevel)
                {
                    // Point levels are different, this can't be a split point,
                    // set flag to false and break immediatelly
                    splitPointCandidate = false;
                    break;
                }
                // else: this is still potential split point candidate so
                //       there's nothing to do
            } // End for all points of this face

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

    // Note: If there is no dynamic load balancing, points at the boundary can't
    // be split points by definition of refinement pattern. However, if there is
    // dynamic load balancing, it may be possible that the split point ends up
    // at the boundary. In that case, the split point should be correctly marked
    // on both sides and unrefined, but this has not been tested. In case of
    // problems, look here first. VV, 26/Jan/2018.

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

    // Extend protected cells across face neighbours
    extendMarkedCells(protectedCell);

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


    // PART 4: Ensure face consistent (2:1 constraint) and possibly point
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

        if (pointBasedRefinement_)
        {
            // Check for 4:1 point based consistent unrefinement. Updates
            // cellsToUnrefine and returns number of removed cells from
            // unrefinement in this iteration
            nRemCells +=
                pointConsistentUnrefinement
                (
                    splitPointsToUnrefine,
                    cellsToUnrefine
                );
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

    Info<< "polyhedralRefinement::setSplitPointsToUnrefine"
        << "(const labelList& unrefinementPointCandidates)" << nl
        << "Selected "
        << returnReduce(splitPointsToUnrefine_.size(), sumOp<label>())
        << " split points to unrefine." << endl;
}


bool Foam::polyhedralRefinement::changeTopology() const
{
    // If either cellsToRefine_ or splitPointsToUnrefine_ have at least one
    // entry, the topology is changing
    if (!cellsToRefine_.empty() || !splitPointsToUnrefine_.empty())
    {
        return true;
    }

    // Otherwise, there is no topo change (at least on this processor in this
    // time step): return false
    return false;
}


void Foam::polyhedralRefinement::setRefinement(polyTopoChange& ref) const
{
    if (!cellsToRefine_.empty())
    {
        // There are some cells to refine, insert the refinement instruction
        setPolyhedralRefinement(ref);
    }

    if (!splitPointsToUnrefine_.empty())
    {
        // There are some split points to unrefine around, insert the
        // unrefinement instruction
        setPolyhedralUnrefinement(ref);
    }

    // Clear the list of cells to refine and split points to unrefine since the
    // refinement/unrefinement instructions have been set
    cellsToRefine_.clear();
    splitPointsToUnrefine_.clear();
}


void Foam::polyhedralRefinement::modifyMotionPoints
(
    pointField& motionPoints
) const
{
    if (debug)
    {
        Pout<< "void polyhedralRefinement::modifyMotionPoints("
            << "pointField& motionPoints) const for object "
            << name() << " : ";
    }

    if (debug)
    {
        Pout << "No motion point adjustment" << endl;
    }
}



void Foam::polyhedralRefinement::updateMesh(const mapPolyMesh& map)
{
    if (debug)
    {
        Info<< "polyhedralRefinement::updateMesh(const mapPolyMesh&) "
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
                    "polyhedralRefinement::updateMesh(const mapPolyMesh& map)"
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

    // Transfer the new point level into the data member
    pointLevel_.transfer(newPointLevel);

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Update face remover
    faceRemover_.updateMesh(map);
}


void Foam::polyhedralRefinement::write(Ostream& os) const
{
    os  << nl << type() << nl
        << name() << nl
        << maxCells_ << nl
        << maxRefinementLevel_ << nl
        << pointBasedRefinement_ << nl
        << nBufferLayers_ << endl;
}


void Foam::polyhedralRefinement::writeDict(Ostream& os) const
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
        << "    pointBasedRefinement " << pointBasedRefinement_
        << token::END_STATEMENT << nl
        << "    nBufferLayers " << nBufferLayers_
        << token::END_STATEMENT << nl
        << "    active " << active()
        << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// ************************************************************************* //
