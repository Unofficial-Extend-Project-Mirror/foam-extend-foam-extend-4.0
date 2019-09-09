/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Hrvoje Jasak, Wikki Ltd.

Notes
    Generalisation of hexRef8 for polyhedral cells and refactorisation into mesh
    modifier engine.

\*---------------------------------------------------------------------------*/

#include "refinement.H"
#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "meshTools.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(refinement, 0);

    // Note: do not add to run-time selection table since this is abstract base
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::refinement::setInstance(const fileName& inst) const
{
    if (debug)
    {
        Pout<< "refinement::setInstance(const fileName& inst)"
            << nl
            << "Resetting file instance of refinement data to " << inst
            << endl;
    }

    cellLevel_.instance() = inst;
    pointLevel_.instance() = inst;
}


Foam::label Foam::refinement::addFace
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


Foam::label Foam::refinement::addInternalFace
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


void Foam::refinement::modifyFace
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


void Foam::refinement::walkFaceToMid
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


void Foam::refinement::walkFaceFromMid
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


Foam::label Foam::refinement::findMinLevel(const labelList& f) const
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


Foam::label Foam::refinement::findMaxLevel(const labelList& f) const
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


Foam::label Foam::refinement::countAnchors
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

void Foam::refinement::adjustRefLevel
(
    label& curNewCellLevel,
    const label oldCellI
)
{
    if (oldCellI == -1)
    {
        // This cell is inflated (does not originate from other cell), set
        // cell level to -1
        curNewCellLevel = -1;
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
            curNewCellLevel = cellLevel_[oldCellI] - 1;
        }
        else if (refStatus == UNCHANGED)
        {
            // Cell hasn't been changed during this refinement, copy old
            // cell level
            curNewCellLevel = cellLevel_[oldCellI];
        }
        else if (refStatus == REFINED)
        {
            // Cell has been refined, increment cell level
            curNewCellLevel = cellLevel_[oldCellI] + 1;
        }
        else
        {
            FatalErrorInFunction
                << "Invalid refinement status detected: "
                << refStatus << nl
                << "Old cell index: " << oldCellI << abort(FatalError);
        }
    }
}


void Foam::refinement::checkInternalOrientation
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
            "refinement::checkInternalOrientation(...)"
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
            "refinement::checkInternalOrientation(...)"
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


void Foam::refinement::checkBoundaryOrientation
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
            "refinement::checkBoundaryOrientation(...)"
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
            "refinement::checkBoundaryOrientation(...)"
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


Foam::label Foam::refinement::faceConsistentRefinement
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
        // Bugfix related to PLB: Check whether owner is already marked for
        // refinement. Will allow 2:1 consistency across certain processor faces
        // where we have a new processor boundary. VV, 23/Jan/2019.
        if
        (
            (neiLevel[i] > curOwnLevel)
         && !cellsToRefine[own]
        )
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


Foam::label Foam::refinement::edgeConsistentRefinement
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
    // with arbitrary polyhedral cells).
    // See faceConsistentRefinement for details. VV, 17/Apr/2018.

    // Return number of added cells
    return nAddCells;
}


Foam::label Foam::refinement::faceConsistentUnrefinement
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
                    "label refinement::faceConsistentUnrefinement"
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
            else
            {
                cellsToUnrefine[own] = false;
                ++nRemCells;
            }
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
                    "label refinement::faceConsistentUnrefinement"
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
            else
            {
                cellsToUnrefine[nei] = false;
                ++nRemCells;
            }
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
        // checks for parallel runs
        if (curOwnLevel < neiLevel[i])
        {
            // Owner level is smaller than neighbour level, we must not
            // unrefine owner

            // Check whether the cell has not been marked for unrefinement
            // Note: this "redundancy" check should not be performed if we are
            // running with dynamic load balancing. If an ordinary face with
            // standard consistency (2:1) becomes a processor face with more
            // stringent consistency (1:1), the refinement still remains valid,
            // even though the 1:1 consistency is not achieved for this time
            // step. Doing further refinement will make sure that this does not
            // exceed at least 2:1 consistency and therefore 2:1 edge
            // consistency as well. Instead of issuing a FatalError, issue a
            // Warning and wrap it into debug
            if (!cellsToUnrefine[own])
            {
                if (debug)
                {
                    WarningIn
                    (
                        "label refinement::faceConsistentUnrefinement"
                        "(boolList& cellsToUnrefine)"
                    )   << "Boundary cell not marked for unrefinement,"
                        << " indicating a previous unnoticed problem with"
                        << " unrefinement."
                        << nl
                        << "Owner: " << own
                        << nl
                        << "Owner level: " << curOwnLevel
                        << ", neighbour level: " << neiLevel[i] << nl
                        << "This is probably because the refinement and "
                        << "unrefinement regions are very close." << nl
                        << "Try increasing nUnrefinementBufferLayers. "
                        << nl
                        << "Another possibility is that you are running "
                        << "with Dynamic Load Balancing, in which case "
                        << "this should be fine."
                        << endl;
                }
            }
            else
            {
                cellsToUnrefine[own] = false;
                ++nRemCells;
            }
        }

        // Note: other possibility (that neighbour level is smaller than owner
        // level) is taken into account on the other side automatically
    }

    // Return number of local cells removed from unrefinement
    return nRemCells;
}


Foam::label Foam::refinement::edgeConsistentUnrefinement
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
                        if (debug)
                        {
                            WarningIn
                            (
                                "label refinement::"
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
                                << endl;
                        }
                    }
                    else
                    {
                        cellsToUnrefine[cellI] = false;
                        ++nRemCells;
                    }
                }
                else if (cellJLevel < cellILevel - 1)
                {
                    // Level of cellJ is smaller than level of cellI - 1, cellJ
                    // must be protected from unrefinement

                    // Check whether the cell has not been marked for
                    // unrefinement
                    if (!cellsToUnrefine[cellJ])
                    {
                        if (debug)
                        {
                            WarningIn
                            (
                                "label refinement::"
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
                                << endl;
                        }
                    }
                    else
                    {
                        cellsToUnrefine[cellJ] = false;
                        ++nRemCells;
                    }
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
    // countless variants when dealing with arbitrary polyhedral cells).
    // See faceConsistentRefinement for details. VV, 3/Apr/2018.

    // Return number of removed cells
    return nRemCells;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinement::refinement
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
        labelField(mesh_.nCells(), 0)
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
        labelField(mesh_.nPoints(), 0)
    ),
    refinementLevelIndicator_(0), // Must be empty before setting refinement
    faceRemover_(mesh_, GREAT),   // Merge boundary faces wherever possible
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
    Info<< "refinement::refinement: " << "Created pointLevel and cellLevel"
        << endl;

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
            "refinement::refinement"
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
            "refinement::refinement"
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
            "refinement::refinement"
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
            "refinement::refinement"
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
            "refinement::refinement"
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
            "refinement::refinement"
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
            "refinement::refinement"
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
            << "Make sure that the number of unrefinement buffer layers is "
            << "at least nRefinementBufferLayers + 2" << nl
            << "in order to avoid problems with edge level inconsistency when "
            << "refinement and unrefinement are performed in same iteration."
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::refinement::~refinement()
{}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

bool Foam::refinement::changeTopology() const
{
    if (!active())
    {
        // Modifier is inactive, skip topo change
        if (debug)
        {
            Pout<< "bool refinement::changeTopology() const"
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


void Foam::refinement::setRefinement(polyTopoChange& ref) const
{
    // Set refinement and unrefinement
    this->setRefinementInstruction(ref);
    this->setUnrefinementInstruction(ref);

    // Clear the list of cells to refine and split points to unrefine since the
    // refinement/unrefinement instructions have been set
    cellsToRefine_.clear();
    splitPointsToUnrefine_.clear();
}


void Foam::refinement::modifyMotionPoints
(
    pointField& motionPoints
) const
{
    if (debug)
    {
        Pout<< "void refinement::modifyMotionPoints("
            << "pointField& motionPoints) const for object "
            << name() << " : ";
    }

    if (debug)
    {
        Pout << "No motion point adjustment" << endl;
    }
}



void Foam::refinement::updateMesh(const mapPolyMesh& map)
{
    if (debug)
    {
        Info<< "refinement::updateMesh(const mapPolyMesh&) "
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
        adjustRefLevel(newCellLevel[newCellI], cellMap[newCellI]);
    }

    // Do cells from points: unrefinement
    // Note: should other types be done as well?
    // HJ, 9/Sep/2019
    const List<objectMap>& cellsFromPoints =  map.cellsFromPointsMap();

    forAll (cellsFromPoints, cfpI)
    {
        adjustRefLevel
        (
            newCellLevel[cellsFromPoints[cfpI].index()],
            cellsFromPoints[cfpI].masterObjects()[0]
        );
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

    // Sync the new point level. Note: here, we assume that the call to
    // updateMesh happened after polyBoundaryMesh::updateMesh where the
    // processor data is fully rebuilt. If this is not the case, the point
    // level will remain unsynced and will cause all kinds of trouble that
    // are extremely difficult to spot. See the change in
    // polyTopoChanger::changeMesh order of calling polyMesh::updateMesh and
    // polyTopoChanger::update. VV, 19/Feb/2019.
    syncTools::syncPointList
    (
        mesh_,
        newPointLevel,
        maxEqOp<label>(),
        label(0),   // Null value
        true        // Apply separation for parallel cyclics
    );

    // Transfer the new point level into the data member
    pointLevel_.transfer(newPointLevel);

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Update face remover
    faceRemover_.updateMesh(map);
}


void Foam::refinement::write(Ostream& os) const
{
    os  << nl << type() << nl
        << name() << nl
        << maxCells_ << nl
        << maxRefinementLevel_ << nl
        << edgeBasedConsistency_ << nl
        << nRefinementBufferLayers_ << nl
        << nUnrefinementBufferLayers_ << endl;
}


void Foam::refinement::writeDict(Ostream& os) const
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
