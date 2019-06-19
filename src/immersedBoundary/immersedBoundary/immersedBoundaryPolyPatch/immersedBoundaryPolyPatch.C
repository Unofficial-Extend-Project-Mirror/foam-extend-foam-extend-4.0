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

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryPolyPatch.H"
#include "foamTime.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "emptyPolyPatch.H"
#include "ImmersedFace.H"
#include "ImmersedCell.H"
#include "triSurfaceDistance.H"
#include "mergePoints.H"
#include "processorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, immersedBoundaryPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        immersedBoundaryPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::debug::tolerancesSwitch
Foam::immersedBoundaryPolyPatch::spanFactor_
(
    "immersedBoundarySpanFactor",
    20
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::immersedBoundaryPolyPatch::cellSpan
(
    const label cellID
) const
{
    const polyMesh& mesh = boundaryMesh().mesh();

    // Calculate span from the bounding box size (prefactor is arbitrary, IG
    // 10/Nov/2018)
    const scalar delta = spanFactor_()*cmptMax
    (
        boundBox
        (
            mesh.cells()[cellID].points
            (
                mesh.faces(),
                mesh.points()
            ),
            false // Do not reduce
        ).span()
    );

    return vector(delta, delta, delta);
}


void Foam::immersedBoundaryPolyPatch::calcTriSurfSearch() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryPolyPatch::calcTriSurfSearch() const")
            << "creating triSurface search algorithm"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already
    if (triSurfSearchPtr_)
    {
        FatalErrorIn
        (
            "void immersedBoundaryPolyPatch::calcTriSurfSearch() const"
        )   << "triSurface search algorithm already exist"
            << abort(FatalError);
    }

    triSurfSearchPtr_ = new triSurfaceSearch(ibMesh_);
}


void Foam::immersedBoundaryPolyPatch::calcImmersedBoundary() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryPolyPatch::calcImmersedBoundary() const")
            << "Calculating geometry"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already
    if
    (
        ibPatchPtr_
     || ibCellsPtr_
     || ibCellCentresPtr_
     || ibCellVolumesPtr_
     || ibFacesPtr_
     || ibFaceCentresPtr_
     || ibFaceAreasPtr_
     || nearestTriPtr_
     || deadCellsPtr_
     || deadFacesPtr_
    )
    {
        FatalErrorIn("immersedBoundaryPolyPatch::calcImmersedBoundary() const")
            << "Geometry already calculated"
            << abort(FatalError);
    }

    // Get reference to the mesh
    const polyBoundaryMesh& bMesh = boundaryMesh();
    const polyMesh& mesh = bMesh.mesh();

    // Get triSurface search
    const triSurfaceSearch& tss = triSurfSearch();

    // Get mesh points
    const pointField& p = mesh.points();

    // Get mesh faces
    const faceList& f = mesh.faces();

    // Get mesh face centres
    const vectorField& Cf = mesh.faceCentres();

    // Get mesh face areas
    const vectorField& S = mesh.faceAreas();

    // Get mesh cell centres
    const vectorField& C = mesh.cellCentres();

    // Get mesh cell volumes
    const scalarField& V = mesh.cellVolumes();

    // Get face addressing
    const labelList& owner = mesh.faceOwner();
    const labelList& neighbour = mesh.faceNeighbour();

    // Get cell-point addressing
    const labelListList& cellPoints = mesh.cellPoints();

    // Algorithm
    // Initialise the search by marking the inside points using calcInside
    // Based on inside points addressing, check intersected faces and cells
    // For all intersected cells, calculate the actual intersection and
    // - calculate the (cell) intersection face, its centre, and area vector
    // - adjust the cell volume and centre
    // - adjust the face area and face centre

    // Mark points that are inside or outside of the triangular surface
    boolList pointsInside = tss.calcInside(p);

    // Adjust selection of points: inside or outside of immersed boundary
    if (internalFlow())
    {
        Info<< "Internal flow" << endl;
    }
    else
    {
        Info<< "External flow" << endl;

        // Flip all points inside identifier
        forAll (pointsInside, i)
        {
            pointsInside[i] = !pointsInside[i];
        }
    }

    // Check cell intersections
    labelList intersectedCell(mesh.nCells(), immersedPoly::UNKNOWN);

    // Estimate the number of intersected cells.
    // Used for sizing of dynamic list only
    // HJ, 11/Dec/2017
    label nIntersectedCells = 0;

    // Go through the faces at the interface between a live and dead cell
    // and mark the band of possible intersections
    forAll (intersectedCell, cellI)
    {
        // Get current cell points
        const labelList& curCp = cellPoints[cellI];

        bool foundInside = false;
        bool foundOutside = false;

        forAll (curCp, cpI)
        {
            if (pointsInside[curCp[cpI]])
            {
                // Found a point inside
                foundInside = true;
            }
            else
            {
                // Found a points outside
                foundOutside = true;
            }
        }

        // Check cell classification
        if (foundInside && !foundOutside)
        {
            // All points inside: cell is wet
            intersectedCell[cellI] = immersedPoly::WET;
        }
        else if (!foundInside && foundOutside)
        {
            // All points outside: cell is dry
            intersectedCell[cellI] = immersedPoly::DRY;
        }
        else if (foundInside && foundOutside)
        {
            // Get span
            const vector span = cellSpan(cellI);

            // If the nearest triangle cannot be found within span than this is
            // most probably a tri surface search error. Mark unknown and check
            // later. (IG 22/Nov/2018)
            if (tss.nearest(C[cellI], span/spanFactor_()).index() == -1)
            {
                intersectedCell[cellI] = immersedPoly::UNKNOWN;
            }
            else
            {
                // Intersected cell
                intersectedCell[cellI] = immersedPoly::CUT;
                nIntersectedCells++;
            }
        }
    }

    // Do a check of the cells selected for cutting but not within the span of
    // the tri surface. The cause of this can either be a stl that is not
    // perfect or ther was an error in the inside/outside tri-search for other
    // reasons. Look at the neigbours that are not CUT and assign their status.
    const cellList& cells = mesh.cells();

    forAll (intersectedCell, cellI)
    {
        if (intersectedCell[cellI] == immersedPoly::UNKNOWN)
        {
            // Check the neigbours
            const cell& curCell = cells[cellI];
            Switch foundWetNei = false;
            Switch foundDryNei = false;

            forAll (curCell, faceI)
            {
                // Only do the check for internal faces. If the face is boundary
                // face then there is nothing to do.
                // NOTE: parallelisation needed?
                if (mesh.isInternalFace(curCell[faceI]))
                {
                    label own = intersectedCell[owner[curCell[faceI]]];
                    label nei = intersectedCell[neighbour[curCell[faceI]]];

                    if
                    (
                        (nei == immersedPoly::DRY)
                     || (own == immersedPoly::DRY)
                    )
                    {
                        foundDryNei = true;
                    }
                    if
                    (
                        (nei == immersedPoly::WET)
                     || (own == immersedPoly::WET)
                    )
                    {
                        foundWetNei = true;
                    }
                }
            }

            if (foundWetNei && !foundDryNei)
            {
                intersectedCell[cellI] = immersedPoly::WET;
            }
            else if (!foundWetNei && foundDryNei)
            {
                intersectedCell[cellI] = immersedPoly::DRY;
            }
            else
            {
                // There are either no wet or dry negbours or there are both.
                // This should not be possible. NOTE: the check is not
                // parallelised and this can theoretically lead to failures in
                // strange arrangaments.
                // Issue a warning, mark CUT and hope for the best.
                // (IG 22/Nov/2018)
                WarningIn
                (
                    "immersedBoundaryPolyPatch::calcImmersedBoundary() const"
                )
                    << "Cannot find wet or dry neigbours! Cell C:"
                    << C[cellI]
                    << " Neighbours: WET:" << foundWetNei
                    << ", DRY:" << foundDryNei
                    << endl;

                intersectedCell[cellI] = immersedPoly::CUT;
                nIntersectedCells++;
            }
        }
    }

    // Count all IB cells and faces for debug
    labelList totalIbCount(4);

    // Collect intersection points and faces.  Primitive patch will be created
    // after renumbering

    // IB points
    // Note: it is difficult to estimate the correct size, so use a guessed
    // number of intersected cells and a dynamic list for automatic resizing
    // HJ, 11/Dec/2017
    DynamicList<point> unmergedPoints
    (
        nIntersectedCells*primitiveMesh::pointsPerFace_
    );
    label nIbPoints = 0;

    // IB patch faces: Cell intersections with the IB patch
    faceList unmergedFaces(mesh.nCells());

    // IB cells: cells intersected by the IB patch
    // This also corresponds to faceCells next to the IB patch
    ibCellsPtr_ = new labelList(mesh.nCells());
    labelList& ibCells = *ibCellsPtr_;

    // IB cellCentres: centre of live part of the intersected cell
    // next to the IB patch
    ibCellCentresPtr_ = new vectorField(mesh.nCells());
    vectorField& ibCellCentres = *ibCellCentresPtr_;

    // IB cellCentres: centre of live part of the intersected cell
    // next to the IB patch
    ibCellVolumesPtr_ = new scalarField(mesh.nCells());
    scalarField& ibCellVolumes = *ibCellVolumesPtr_;

    // Nearest triangle
    nearestTriPtr_ = new labelList(mesh.nCells());
    labelList& nearestTri = *nearestTriPtr_;

    // Count interected cells
    label nIbCells = 0;

    // Collect dead cells
    boolList deadCells(mesh.nCells(), false);

    // At this point, all live cells are marked with 1
    // Intesect all cells that are marked for intersection

    forAll (intersectedCell, cellI)
    {
        if (intersectedCell[cellI] == immersedPoly::CUT)
        {
            // Found intersected cell

            // Get span
            const vector span = cellSpan(cellI);

            // Create a cutting object with a local tolerance
            triSurfaceDistance dist
            (
                tss,
                2*span,
                internalFlow(),
                true               // iterate intersection
            );

            // Calculate the intersection
            ImmersedCell<triSurfaceDistance> cutCell
            (
                cellI,
                mesh,
                dist
            );

            // Check for irregular intersections
            if (cutCell.isAllWet())
            {
                intersectedCell[cellI] = immersedPoly::WET;
            }
            else if (cutCell.isAllDry())
            {
                intersectedCell[cellI] = immersedPoly::DRY;
            }
            else
            {
                // True intersection.  Cut the cell and store all
                // derived data

                // Note: volumetric check is not allowed because true
                // intersection guarantees that the faces of the cell
                // have been cut.  Therefore, the cell MUST be an IB cell.
                // If the cut is invalid, Marooney Maneouvre shall correct
                // the error in sum(Sf).  HJ, 12/Mar/2019

                // Store ibFace with local points. Points merge will
                // take place later
                const face& cutFace = cutCell.faces()[0];

                const pointField& cutPoints = cutCell.points();

                // Collect the renumbered face, using the point labels
                // from the unmergedPoints list
                face renumberedFace(cutFace.size());

                // Insert points and renumber the face
                forAll (cutFace, cpI)
                {
                    unmergedPoints.append(cutPoints[cutFace[cpI]]);
                    renumberedFace[cpI] = nIbPoints;
                    nIbPoints++;
                }

                // Record the face
                unmergedFaces[nIbCells] = renumberedFace;

                // Collect cut cell index
                ibCells[nIbCells] = cellI;

                // Record the live centre
                ibCellCentres[nIbCells] = cutCell.wetVolumeCentre();

                // Record the live volume
                ibCellVolumes[nIbCells] = cutCell.wetVolume();

                // Record the nearest triangle to the face centre
                nearestTri[nIbCells] =
                    tss.nearest(cutFace.centre(cutPoints), span).index();

                nIbCells++;
            }
        }
    }

    // Pick up direct face cuts after regular cell cuts are collected
    forAll (neighbour, faceI)
    {
        if
        (
            intersectedCell[owner[faceI]] == immersedPoly::WET
         && intersectedCell[neighbour[faceI]] == immersedPoly::DRY
        )
        {
            // Direct face cut, owner

            // Grab a point and wet cell and make an IB face
            pointField facePoints = f[faceI].points(p);
            face renumberedFace(facePoints.size());

            // Insert points
            forAll (facePoints, fpI)
            {
                unmergedPoints.append(facePoints[fpI]);
                renumberedFace[fpI] = nIbPoints;
                nIbPoints++;
            }

            // Record the face
            unmergedFaces[nIbCells] = renumberedFace;

            // Collect cut cell index
            ibCells[nIbCells] = owner[faceI];

            // Record the live centre
            ibCellCentres[nIbCells] = C[owner[faceI]];

            // Record the live volume: equal to owner volume
            ibCellVolumes[nIbCells] = V[owner[faceI]];

            // Get span of owner and neighbour
            vector span = cellSpan(owner[faceI]);

            span = Foam::max
            (
                span,
                cellSpan(neighbour[faceI])
            );

            // Record the nearest triangle to the face centre
            nearestTri[nIbCells] = tss.nearest(Cf[faceI], span).index();

            nIbCells++;
        }
        else if
        (
            intersectedCell[owner[faceI]] == immersedPoly::DRY
         && intersectedCell[neighbour[faceI]] == immersedPoly::WET
        )
        {
            // Direct face cut, neighbour

            // Grab a point and wet cell and make an IB face
            // Note: reverse face in cut
            pointField facePoints = f[faceI].reverseFace().points(p);

            face renumberedFace(facePoints.size());

            // Insert points
            forAll (facePoints, fpI)
            {
                unmergedPoints.append(facePoints[fpI]);
                renumberedFace[fpI] = nIbPoints;
                nIbPoints++;
            }

            // Record the face
            unmergedFaces[nIbCells] = renumberedFace;

            // Collect cut cell index
            ibCells[nIbCells] = neighbour[faceI];

            // Record the live centre
            ibCellCentres[nIbCells] = C[neighbour[faceI]];

            // Record the live volume: equal to neighbour volume
            ibCellVolumes[nIbCells] = V[neighbour[faceI]];

            // Get span of neighbour and neighbour
            vector span = cellSpan(neighbour[faceI]);

            span = Foam::max
            (
                span,
                cellSpan(owner[faceI])
            );

            // Record the nearest triangle to the face centre
            nearestTri[nIbCells] = tss.nearest(Cf[faceI], span).index();

            nIbCells++;
        }
    }

    // Check coupled boundaries for direct face cuts

    // Assemble local and neighbour cuts for coupled patches only
    labelListList coupledPatchOwnCut(bMesh.size());
    labelListList coupledPatchNbrCut(bMesh.size());

    // Note: this part requires a rewrite using virtual functions
    // to communicate the cut data from the shadow cell
    // (across the coupled interface) in order to determine
    // the coupled face status.
    // Currently, this is enabled only for processor boundaries.
    // HJ, 28/Dec/2017

    // Send loop
    forAll (bMesh, patchI)
    {
        if (bMesh[patchI].coupled())
        {
            if (isA<processorPolyPatch>(bMesh[patchI]))
            {
                if (Pstream::parRun())
                {
                    const processorPolyPatch& curProcPatch =
                        refCast<const processorPolyPatch>(bMesh[patchI]);

                    // Send internal cut
                    coupledPatchOwnCut[patchI] = labelList
                    (
                        intersectedCell,
                        bMesh[patchI].faceCells()
                    );

                    OPstream toNeighbProc
                    (
                        Pstream::blocking,
                        curProcPatch.neighbProcNo(),
                        sizeof(label)*curProcPatch.size()
                    );

                    toNeighbProc << coupledPatchOwnCut[patchI];
                }
            }
            else
            {
                WarningIn
                (
                    "void immersedBoundaryPolyPatch::"
                    "calcImmersedBoundary() const"
                )   << "Non-processor coupled patch detected for "
                    << "immersed boundary.  "
                    << "Direct face cut may not be detected"
                    << endl;
            }
        }
    }

    // Receive loop
    forAll (bMesh, patchI)
    {
        if (bMesh[patchI].coupled())
        {
            if (isA<processorPolyPatch>(bMesh[patchI]))
            {
                if (Pstream::parRun())
                {
                    const processorPolyPatch& curProcPatch =
                        refCast<const processorPolyPatch>(bMesh[patchI]);

                    IPstream fromNeighbProc
                    (
                        Pstream::blocking,
                        curProcPatch.neighbProcNo(),
                        sizeof(label)*curProcPatch.size()
                    );

                    coupledPatchNbrCut[patchI] = labelList(fromNeighbProc);
                }
            }
        }
    }

    // Analyse the cut
    forAll (bMesh, patchI)
    {
        if (!coupledPatchOwnCut[patchI].empty())
        {
            const labelList& curOwnCut = coupledPatchOwnCut[patchI];
            const labelList& curNbrCut = coupledPatchNbrCut[patchI];

            const labelList& fc = bMesh[patchI].faceCells();

            forAll (curOwnCut, patchFaceI)
            {
                if
                (
                    curOwnCut[patchFaceI] == immersedPoly::WET
                 && curNbrCut[patchFaceI] == immersedPoly::DRY
                )
                {
                    // Direct face cut, coupled on live side

                    // Get face index.  Note the difference between faceI
                    // and patchFaceI
                    const label faceI = bMesh[patchI].start() + patchFaceI;

                    // Grab a point and wet cell and make an IB face
                    pointField facePoints = f[faceI].points(p);
                    face renumberedFace(facePoints.size());

                    // Insert points
                    forAll (facePoints, fpI)
                    {
                        unmergedPoints.append(facePoints[fpI]);
                        renumberedFace[fpI] = nIbPoints;
                        nIbPoints++;
                    }

                    // Record the face
                    unmergedFaces[nIbCells] = renumberedFace;

                    // Collect cut cell index
                    ibCells[nIbCells] = fc[patchFaceI];

                    // Record the live centre
                    ibCellCentres[nIbCells] = C[fc[patchFaceI]];

                    // Record the live volume: equal to owner volume
                    ibCellVolumes[nIbCells] = V[fc[patchFaceI]];

                    // Get span of owner.  Cannot reach neighbour
                    vector span = cellSpan(fc[patchFaceI]);

                    // Record the nearest triangle to the face centre
                    nearestTri[nIbCells] =
                        tss.nearest(Cf[faceI], span).index();

                    nIbCells++;
                }
            }
        }
    }

    // Record the number of IB cells for debug
    totalIbCount[0] = nIbCells;

    // Reset the cell lists
    unmergedFaces.setSize(nIbCells);
    ibCells.setSize(nIbCells);
    ibCellCentres.setSize(nIbCells);
    ibCellVolumes.setSize(nIbCells);
    nearestTri.setSize(nIbCells);

    // Check tri addressing
    if (min(nearestTri) == -1)
    {
        FatalErrorInFunction
            << "Cannot find nearestTri for all points"
            << abort(FatalError);
    }

    // Build stand-alone patch
    // Memory management
    {
        unmergedPoints.shrink();

        pointField ibPatchPoints;
        labelList pointMap;

        mergePoints
        (
            unmergedPoints,
            1e-6,               // mergeTol.  Review.  Do not like the algorithm
            false,              // verbose
            pointMap,
            ibPatchPoints
        );

        // Renumber faces after point merge
        faceList ibPatchFaces(unmergedFaces.size());

        forAll (unmergedFaces, faceI)
        {
            // Get old and new face
            const face& uFace = unmergedFaces[faceI];
            face& rFace = ibPatchFaces[faceI];
            rFace.setSize(uFace.size());
            forAll (uFace, pointI)
            {
                rFace[pointI] = pointMap[uFace[pointI]];
            }
        }

        // Create IB patch from renumbered points and faces
        ibPatchPtr_ = new standAlonePatch(ibPatchFaces, ibPatchPoints);

        if (mesh.time().outputTime())
        {
            Info << "Writing immersed patch as VTK" << endl;

            fileName fvPath(mesh.time().path()/"VTK");
            mkDir(fvPath);

            fileName surfaceFileName
            (
                "immersed" + name() + "_live_"
              + Foam::name(boundaryMesh().mesh().time().timeIndex())
            );

            ibPatchPtr_->writeVTK(fvPath/surfaceFileName);

            fileName normalsFileName
            (
                "normals" + name() + "_live_"
              + Foam::name(boundaryMesh().mesh().time().timeIndex())
            );

            ibPatchPtr_->writeVTKNormals(fvPath/normalsFileName);
        }
    }

    // Count and collect dead cells

    // Memory management
    {
        label nDeadCells = 0;

        forAll (intersectedCell, cellI)
        {
            if (intersectedCell[cellI] == immersedPoly::DRY)
            {
                nDeadCells++;
            }
        }

        // Allocate storage and collect dead cells
        deadCellsPtr_ = new labelList(nDeadCells);
        labelList& dc = *deadCellsPtr_;

        // Reset the counter
        nDeadCells = 0;

        forAll (intersectedCell, cellI)
        {
            if (intersectedCell[cellI] == immersedPoly::DRY)
            {
                dc[nDeadCells] = cellI;
                nDeadCells++;
            }
        }

        // Record the number of dead cells for debug
        totalIbCount[1] = nDeadCells;
    }

    // IB faces: faces intersected by the IB patch
    // This also corresponds to faceCells next to the IB patch
    ibFacesPtr_ = new labelList(mesh.nFaces());
    labelList& ibFaces = *ibFacesPtr_;

    // IB face centres: centre of live part of the intersected face
    // next to the IB patch
    ibFaceCentresPtr_ = new vectorField(mesh.nFaces());
    vectorField& ibFaceCentres = *ibFaceCentresPtr_;

    // IB face areas: surface-normal area of live part of the intersected face
    // next to the IB patch
    ibFaceAreasPtr_ = new vectorField(mesh.nFaces());
    vectorField& ibFaceAreas = *ibFaceAreasPtr_;
    label nIbFaces = 0;

    // Classify faces
    labelList intersectedFace(mesh.nFaces(), immersedPoly::UNKNOWN);

    // Resolve simple face intersections based on the cell intersection data
    // First, kill all faces touching dead cells, including internal
    // and boundary faces.
    // If a face touches a live cell, it is live
    // The intersection belt will be handled separately by detailed intersection

    // Quick intersection scan: if owner and neighbour are in the same state
    // the face is in the same state

    // Internal faces
    forAll (neighbour, faceI)
    {
        // Wet on wet
        if
        (
            intersectedCell[owner[faceI]] == immersedPoly::WET
         && intersectedCell[neighbour[faceI]] == immersedPoly::WET
        )
        {
            intersectedFace[faceI] = immersedPoly::WET;
        }

        // Dry on dry
        if
        (
            intersectedCell[owner[faceI]] == immersedPoly::DRY
         && intersectedCell[neighbour[faceI]] == immersedPoly::DRY
        )
        {
            intersectedFace[faceI] = immersedPoly::DRY;
        }

        // Wet on cut face must remain wet.  Error in cut cell is fixed
        // by the Marooney Maneouvre.  HJ, 5/Apr/2019
        if
        (
            (
                intersectedCell[owner[faceI]] == immersedPoly::WET
             && intersectedCell[neighbour[faceI]] == immersedPoly::CUT
            )
         || (
                intersectedCell[owner[faceI]] == immersedPoly::CUT
             && intersectedCell[neighbour[faceI]] == immersedPoly::WET
            )
        )
        {
            intersectedFace[faceI] = immersedPoly::WET;
        }

        // Special check for directly cut faces
        // Wet-to-dry and dry-to-wet is a direct face cut
        // Dry-to-cut or cut-to-dry are cutting errors.  They will be
        // corrected later in corrected face areas, based on closed cell
        // tolerance.  HJ, 11/Dec/2017
        if
        (
            (
                intersectedCell[owner[faceI]] == immersedPoly::WET
             && intersectedCell[neighbour[faceI]] == immersedPoly::DRY
            )
         || (
                intersectedCell[owner[faceI]] == immersedPoly::DRY
             && intersectedCell[neighbour[faceI]] == immersedPoly::WET
             )
         || (
                intersectedCell[owner[faceI]] == immersedPoly::DRY
             && intersectedCell[neighbour[faceI]] == immersedPoly::CUT
             )
         || (
                intersectedCell[owner[faceI]] == immersedPoly::CUT
             && intersectedCell[neighbour[faceI]] == immersedPoly::DRY
             )
        )
        {
            // Note:
            // Wet-to-dry: this face has been declared to be a
            //             cut face and needs to be taken out as live face
            // Cut-to-dry: this is either an outside edge of cut faces or
            //             a cutting error
            intersectedFace[faceI] = immersedPoly::DRY;
        }
    }

    // Boundary faces
    forAll (bMesh, patchI)
    {
        const label patchStart = bMesh[patchI].start();

        if (bMesh[patchI].coupled())
        {
            // Coupled patch: two-sided check
            const labelList& curOwnCut = coupledPatchOwnCut[patchI];
            const labelList& curNbrCut = coupledPatchNbrCut[patchI];

            forAll (curOwnCut, patchFaceI)
            {
                // Wet on wet
                if
                (
                    curOwnCut[patchFaceI] == immersedPoly::WET
                 && curNbrCut[patchFaceI] == immersedPoly::WET
                )
                {
                    intersectedFace[patchStart + patchFaceI] =
                        immersedPoly::WET;
                }

                // Dry on dry
                if
                (
                    curOwnCut[patchFaceI] == immersedPoly::DRY
                 && curNbrCut[patchFaceI] == immersedPoly::DRY
                )
                {
                    intersectedFace[patchStart + patchFaceI] =
                        immersedPoly::DRY;
                }

                // Wet on cut face must remain wet.  Error in cut cell is fixed
                // by the Marooney Maneouvre.  HJ, 5/Apr/2019
                if
                (
                    (
                        curOwnCut[patchFaceI] == immersedPoly::WET
                     && curNbrCut[patchFaceI] == immersedPoly::CUT
                    )
                 || (
                        curOwnCut[patchFaceI] == immersedPoly::CUT
                     && curNbrCut[patchFaceI] == immersedPoly::WET
                    )
                )
                {
                    intersectedFace[patchStart + patchFaceI] =
                        immersedPoly::WET;
                }

                // Special check for directly cut faces
                // Wet-to-dry and dry-to-wet is a direct face cut
                // Dry-to-cut or cut-to-dry are cutting errors.  They will be
                // corrected later in corrected face areas, based on closed cell
                // tolerance.  HJ, 11/Dec/2017
                if
                (
                    (
                        curOwnCut[patchFaceI] == immersedPoly::WET
                      && curNbrCut[patchFaceI] == immersedPoly::DRY
                    )
                 || (
                        curOwnCut[patchFaceI] == immersedPoly::DRY
                     && curNbrCut[patchFaceI] == immersedPoly::WET
                    )
                 || (
                        curOwnCut[patchFaceI] == immersedPoly::DRY
                     && curNbrCut[patchFaceI] == immersedPoly::CUT
                    )
                 || (
                        curOwnCut[patchFaceI] == immersedPoly::CUT
                     && curNbrCut[patchFaceI] == immersedPoly::DRY
                     )
                )
                {
                    // Note:
                    // Wet-to-dry: this face has been declared to be a
                    // cut face and needs to be taken out as live face
                    // Cut-to-dry: this is either an outside edge of cut faces
                    // or a cutting error
                    intersectedFace[patchStart + patchFaceI] =
                        immersedPoly::DRY;
                }
            }
        }
        else
        {
            // Regular patch: one-sided check
            const labelList& fc = bMesh[patchI].faceCells();

            forAll (fc, patchFaceI)
            {
                if
                (
                    intersectedCell[fc[patchFaceI]] == immersedPoly::WET
                )
                {
                    intersectedFace[patchStart + patchFaceI] =
                        immersedPoly::WET;
                }

                if
                (
                    intersectedCell[fc[patchFaceI]] == immersedPoly::DRY
                )
                {
                    intersectedFace[patchStart + patchFaceI] =
                        immersedPoly::DRY;
                }
            }
        }
    }

    // Detailed face check after initial rejection scan
    forAll (intersectedFace, faceI)
    {
        if (intersectedFace[faceI] == immersedPoly::UNKNOWN)
        {
            // Possibly intersected face.  Check existance of intersection
            // via points
            const labelList& curF = f[faceI];

            bool foundInside = false;
            bool foundOutside = false;

            forAll (curF, fI)
            {
                if (pointsInside[curF[fI]])
                {
                    // Found a point inside
                    foundInside = true;
                }
                else
                {
                    // Found a points outside
                    foundOutside = true;
                }
            }

            // Check face classification
            if (foundInside && !foundOutside)
            {
                // All points inside: cell is wet
                intersectedFace[faceI] = immersedPoly::WET;
            }
            else if (!foundInside && foundOutside)
            {
                // All points outside: cell is dry
                intersectedFace[faceI] = immersedPoly::DRY;
            }
            else if (foundInside && foundOutside)
            {
                // Real intersection.  Try to cut the face

                // Get search span
                vector span = cellSpan(owner[faceI]);

                // For internal face, check the neighbour span as well
                if (mesh.isInternalFace(faceI))
                {
                    span = Foam::max
                    (
                        span,
                        cellSpan(neighbour[faceI])
                    );
                }

                // Create a cutting object with a local tolerance
                triSurfaceDistance dist
                (
                    tss,
                    span,
                    internalFlow(),
                    true               // iterate intersection
                );

                // Calculate the intersection
                ImmersedFace<triSurfaceDistance> cutFace
                (
                    faceI,
                    mesh,
                    dist
                );

                if (cutFace.isAllWet())
                {
                    intersectedFace[faceI] = immersedPoly::WET;
                }
                else if (cutFace.isAllDry())
                {
                    intersectedFace[faceI] = immersedPoly::DRY;
                }
                else
                {
                    // Real intesection.  Check cut. Rejection on thin cut is
                    // performed by ImmersedFace.  HJ, 13/Mar/2019
                    const scalar faceFactor =
                        cutFace.wetAreaMag()/mag(S[faceI]);

                    // True intersection.  Collect data
                    intersectedFace[faceI] = immersedPoly::CUT;

                    // Get intersected face index
                    ibFaces[nIbFaces] = faceI;

                    // Get wet centre
                    ibFaceCentres[nIbFaces] = cutFace.wetAreaCentre();

                    // Get wet area, preserving original normal direction
                    ibFaceAreas[nIbFaces] = faceFactor*S[faceI];

                    nIbFaces++;
                }
            }
        }
    }

    // Record the number of IB faces for debug
    totalIbCount[2] = nIbFaces;

    // Reset the sizes of the list
    ibFaces.setSize(nIbFaces);
    ibFaceCentres.setSize(nIbFaces);

    // Count and collect dead faces
    // Memory management
    {
        label nDeadFaces = 0;

        forAll (intersectedFace, faceI)
        {
            if (intersectedFace[faceI] == immersedPoly::DRY)
            {
                nDeadFaces++;
            }
        }

        // Allocate storage and collect dead faces
        deadFacesPtr_ = new labelList(nDeadFaces);
        labelList& df = *deadFacesPtr_;

        // Reset the counter
        nDeadFaces = 0;

        forAll (intersectedFace, faceI)
        {
            if (intersectedFace[faceI] == immersedPoly::DRY)
            {
                df[nDeadFaces] = faceI;
                nDeadFaces++;
            }
        }

        // Record the number of dead faces for debug
        totalIbCount[3] = nDeadFaces;
    }

    // Reduce is not allowed in parallel load balancing
    // HJ, 24/Oct/2018
    // if (debug)
    {
        // reduce(totalIbCount, sumOp<List<label> >());

        Info<< "Immersed boundary " << name() << " info: "
            << "nIbCells: " << totalIbCount[0]
            << " nDeadCells: " << totalIbCount[1]
            << " nIbFaces: " << totalIbCount[2]
            << " nDeadFaces: " << totalIbCount[3]
            << endl;
    }
}


void Foam::immersedBoundaryPolyPatch::calcCorrectedGeometry() const
{
    if (debug)
    {
        InfoIn
        (
            "void immersedBoundaryPolyPatch::calcCorrectedGeometry() const"
        )   << "Calculating corrected geometry"
            << endl;
    }

    if
    (
        correctedIbPatchFaceAreasPtr_
    )
    {
        FatalErrorIn
        (
            "void immersedBoundaryPolyPatch::calcCorrectedGeometry() const"
        )   << "Corrected geometry already calculated"
            << abort(FatalError);
    }

    // Get mesh reference
    const polyMesh& mesh = boundaryMesh().mesh();

    // Get mesh geometry references from the MeshObject
    vectorField& C = correctedFields_.correctedCellCentres();;
    vectorField& Cf = correctedFields_.correctedFaceCentres();
    scalarField& V = correctedFields_.correctedCellVolumes();
    vectorField& Sf = correctedFields_.correctedFaceAreas();

    // Initialise IB patch face areas with the areas of the stand-alone patch
    // They will be corrected using the Marooney Maneouvre
    correctedIbPatchFaceAreasPtr_ = new vectorField(ibPatch().areas());
    vectorField& ibSf = *correctedIbPatchFaceAreasPtr_;

    // Correct for all cut cells

    // Get cut cells
    const labelList& cutCells = ibCells();
    const vectorField& cutCellCentres = ibCellCentres();
    const scalarField& cutCellVolumes = ibCellVolumes();

    forAll (cutCells, ccI)
    {
        // Correct the volume and area
        C[cutCells[ccI]] = cutCellCentres[ccI];

        V[cutCells[ccI]] = cutCellVolumes[ccI];
    }

    // Deactivate dead cells
    const labelList& dc = deadCells();

    forAll (dc, dcI)
    {
        // Scale dead volume to small
        V[dc[dcI]] *= SMALL;
    }

    // Correct for all cut faces

    // Get cut faces
    const labelList& cutFaces = ibFaces();
    const vectorField& cutFaceCentres = ibFaceCentres();
    const vectorField& cutFaceAreas = ibFaceAreas();

    forAll (cutFaces, cfI)
    {
        Cf[cutFaces[cfI]] = cutFaceCentres[cfI];

        // Preserve the original face normal
        Sf[cutFaces[cfI]] = cutFaceAreas[cfI];
    }

    // Deactivate dead faces
    const labelList& df = deadFaces();

    forAll (df, dfI)
    {
        // Scale dead area to small
        Sf[df[dfI]] *= SMALL;
    }

    // In case of cutting errors due to finite tolerance, some cut cells may
    // remain opened and have to be closed by force.  This will be achieved
    // by the Marooney Maneouvre, where the face sum imbalance is compensated
    // in the cut face.  HJ, 11/Dec/2017

    const labelList& owner = mesh.faceOwner();

    label nMarooneyCells = 0;

    // Get valid directions to avoid round-off errors in 2-D cases
    const Vector<label> dirs = mesh.geometricD();
    vector validDirs = vector::zero;

    for (direction cmpt = 0; cmpt < Vector<label>::nComponents; cmpt++)
    {
        if (dirs[cmpt] > 0)
        {
            validDirs[cmpt] = 1;
        }
    }

    forAll (cutCells, cutCellI)
    {
        const label ccc = cutCells[cutCellI];

        // Calculate sum Sf and sumMagSf for the cell
        const cell& curCell = mesh.cells()[ccc];

        vector curSumSf = vector::zero;
        scalar curSumMagSf = 0;

        // Collect from regular faces
        forAll (curCell, cfI)
        {
            const vector& curSf = Sf[curCell[cfI]];

            // Check owner/neighbour
            if (owner[curCell[cfI]] == ccc)
            {
                curSumSf += curSf;
            }
            else
            {
                curSumSf -= curSf;
            }

            curSumMagSf += mag(curSf);
        }

        // Add cut face only into mag.  The second part is handled in the
        // if-statement
        curSumMagSf += mag(ibSf[cutCellI]);

        // Adjustment is peformed when the openness is greater than a certain
        // fraction of surface area.  Criterion by IG, 13/Mar/2019
        // Switched to using absolute check from primitiveMeshCheck.
        // HJ, 13/Mar/2019
        // if (mag(curSumSf + ibSf[cutCellI]) > 1e-6*curSumMagSf)
        if (mag(curSumSf + ibSf[cutCellI]) > primitiveMesh::closedThreshold_)
        {
            // Info<< "Marooney Maneouvre for cell " << ccc
            //     << " error: " << curSumSf + ibSf[cutCellI] << " "
            //     << " V: " << cutCellVolumes[cutCellI]
            //     << " Sf: " << ibSf[cutCellI]
            //     << " corr S: " << curSumSf << endl;

            nMarooneyCells++;

            // Create IB face to ideally close the cell
            ibSf[cutCellI] = cmptMultiply(validDirs, -curSumSf);
        }
    }

    if (debug)
    {
        if (nMarooneyCells > 0)
        {
            InfoIn
            (
                "void immersedBoundaryPolyPatch::calcCorrectedGeometry() const"
            )   << "Marooney Maneouvre used for " << nMarooneyCells
                << " out of " << cutCells.size()
                << endl;
        }
    }

    if (min(mag(ibSf)) < SMALL)
    {
        WarningInFunction
            << "Minimum IB face area for patch " << name()
            << ": " << min(mag(ibSf)) << ".  Possible cutting error.  "
            << "Review immersed boundary tolerances."
            << endl;
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::immersedBoundaryPolyPatch::movePoints(const pointField& p)
{
    if (debug)
    {
        InfoIn
        (
            "void immersedBoundaryPolyPatch::"
            "movePoints(const pointField&) const"
        )   << "Moving mesh: immersedBoundary update"
            << endl;
    }

    // Handle motion of the mesh for new immersed boundary position
    if (ibUpdateTimeIndex_ < boundaryMesh().mesh().time().timeIndex())
    {
        // New motion in the current time step.  Clear
        ibUpdateTimeIndex_ = boundaryMesh().mesh().time().timeIndex();

        clearOut();
    }

    polyPatch::movePoints(p);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryPolyPatch::immersedBoundaryPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm),
    ibMesh_
    (
        IOobject
        (
            name  + ".ftr",
            bm.mesh().time().constant(), // instance
            "triSurface",                // local
            bm.mesh().parent(),          // registry
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    internalFlow_(false),
    isWall_(true),
    movingIb_(false),
    ibUpdateTimeIndex_(-1),
    triSurfSearchPtr_(nullptr),
    ibPatchPtr_(nullptr),
    ibCellsPtr_(nullptr),
    ibCellCentresPtr_(nullptr),
    ibCellVolumesPtr_(nullptr),
    ibFacesPtr_(nullptr),
    ibFaceCentresPtr_(nullptr),
    ibFaceAreasPtr_(nullptr),
    nearestTriPtr_(nullptr),
    deadCellsPtr_(nullptr),
    deadFacesPtr_(nullptr),
    correctedFields_(immersedBoundaryCorrectedMeshFields::New(bm.mesh())),
    correctedIbPatchFaceAreasPtr_(nullptr),
    topoChangeIndex_(-1),
    oldIbPointsPtr_(nullptr)
{}


Foam::immersedBoundaryPolyPatch::immersedBoundaryPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, dict, index, bm),
    ibMesh_
    (
        IOobject
        (
            name  + ".ftr",
            bm.mesh().time().constant(), // instance
            "triSurface",                // local
            bm.mesh().parent(),          // read from parent registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    internalFlow_(dict.lookup("internalFlow")),
    isWall_(dict.lookup("isWall")),
    movingIb_(false),
    ibUpdateTimeIndex_(-1),
    triSurfSearchPtr_(nullptr),
    ibPatchPtr_(nullptr),
    ibCellsPtr_(nullptr),
    ibCellCentresPtr_(nullptr),
    ibCellVolumesPtr_(nullptr),
    ibFacesPtr_(nullptr),
    ibFaceCentresPtr_(nullptr),
    ibFaceAreasPtr_(nullptr),
    nearestTriPtr_(nullptr),
    deadCellsPtr_(nullptr),
    deadFacesPtr_(nullptr),
    correctedFields_(immersedBoundaryCorrectedMeshFields::New(bm.mesh())),
    correctedIbPatchFaceAreasPtr_(nullptr),
    topoChangeIndex_(-1),
    oldIbPointsPtr_(nullptr)
{
    if (size() > 0)
    {
        FatalIOErrorIn
        (
            "immersedBoundaryPolyPatch::immersedBoundaryPolyPatch\n"
            "(\n"
            "    const word& name,\n"
            "    const dictionary& dict,\n"
            "    const label index,\n"
            "    const polyBoundaryMesh& bm\n"
            ")",
            dict
        )   << "Faces detected in the immersedBoundaryPolyPatch.  "
            << "This is not allowed: please make sure that the patch size "
            << "equals zero."
            << abort(FatalIOError);
    }
}


Foam::immersedBoundaryPolyPatch::immersedBoundaryPolyPatch
(
    const immersedBoundaryPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    ibMesh_
    (
        IOobject
        (
            pp.name() + ".ftr",
            bm.mesh().time().constant(), // instance
            "triSurface",                // local
            bm.mesh().parent(),          // parent registry
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pp.ibMesh()   // Take ibMesh from pp
    ),
    internalFlow_(pp.internalFlow_),
    isWall_(pp.isWall_),
    movingIb_(false),
    ibUpdateTimeIndex_(-1),
    triSurfSearchPtr_(nullptr),
    ibPatchPtr_(nullptr),
    ibCellsPtr_(nullptr),
    ibCellCentresPtr_(nullptr),
    ibCellVolumesPtr_(nullptr),
    ibFacesPtr_(nullptr),
    ibFaceCentresPtr_(nullptr),
    ibFaceAreasPtr_(nullptr),
    nearestTriPtr_(nullptr),
    deadCellsPtr_(nullptr),
    deadFacesPtr_(nullptr),
    correctedFields_(immersedBoundaryCorrectedMeshFields::New(bm.mesh())),
    correctedIbPatchFaceAreasPtr_(nullptr),
    topoChangeIndex_(-1),
    oldIbPointsPtr_(nullptr)
{}


Foam::immersedBoundaryPolyPatch::immersedBoundaryPolyPatch
(
    const immersedBoundaryPolyPatch& pp
)
:
    polyPatch(pp),
    ibMesh_
    (
        IOobject
        (
            pp.name() + ".ftr",
            pp.boundaryMesh().mesh().time().constant(), // instance
            "triSurface",                               // local
            pp.boundaryMesh().mesh().parent(),          // parent registry
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pp.ibMesh()   // Take ibMesh from pp
    ),
    internalFlow_(pp.internalFlow_),
    isWall_(pp.isWall_),
    movingIb_(false),
    ibUpdateTimeIndex_(-1),
    triSurfSearchPtr_(nullptr),
    ibPatchPtr_(nullptr),
    ibCellsPtr_(nullptr),
    ibCellCentresPtr_(nullptr),
    ibCellVolumesPtr_(nullptr),
    ibFacesPtr_(nullptr),
    ibFaceCentresPtr_(nullptr),
    ibFaceAreasPtr_(nullptr),
    nearestTriPtr_(nullptr),
    deadCellsPtr_(nullptr),
    deadFacesPtr_(nullptr),
    correctedFields_
    (
        immersedBoundaryCorrectedMeshFields::New
        (
            pp.boundaryMesh().mesh()
        )
    ),
    correctedIbPatchFaceAreasPtr_(nullptr),
    topoChangeIndex_(-1),
    oldIbPointsPtr_(nullptr)
{}


Foam::immersedBoundaryPolyPatch::immersedBoundaryPolyPatch
(
    const immersedBoundaryPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    ibMesh_
    (
        IOobject
        (
            pp.name() + ".ftr",
            bm.mesh().time().constant(), // instance
            "triSurface",                // local
            bm.mesh().parent(),          // parent registry
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pp.ibMesh()   // Take ibMesh from pp
    ),
    internalFlow_(pp.internalFlow_),
    isWall_(pp.isWall_),
    movingIb_(false),
    ibUpdateTimeIndex_(-1),
    triSurfSearchPtr_(nullptr),
    ibPatchPtr_(nullptr),
    ibCellsPtr_(nullptr),
    ibCellCentresPtr_(nullptr),
    ibCellVolumesPtr_(nullptr),
    ibFacesPtr_(nullptr),
    ibFaceCentresPtr_(nullptr),
    ibFaceAreasPtr_(nullptr),
    nearestTriPtr_(nullptr),
    deadCellsPtr_(nullptr),
    deadFacesPtr_(nullptr),
    correctedFields_(immersedBoundaryCorrectedMeshFields::New(bm.mesh())),
    correctedIbPatchFaceAreasPtr_(nullptr),
    topoChangeIndex_(-1),
    oldIbPointsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::immersedBoundaryPolyPatch::~immersedBoundaryPolyPatch()
{
    clearOut();

    deleteDemandDrivenData(oldIbPointsPtr_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::triSurfaceSearch&
Foam::immersedBoundaryPolyPatch::triSurfSearch() const
{
    if (!triSurfSearchPtr_)
    {
        calcTriSurfSearch();
    }

    return *triSurfSearchPtr_;
}

const Foam::standAlonePatch&
Foam::immersedBoundaryPolyPatch::ibPatch() const
{
    if (!ibPatchPtr_)
    {
        calcImmersedBoundary();
    }

    return *ibPatchPtr_;
}


const Foam::labelList&
Foam::immersedBoundaryPolyPatch::ibCells() const
{
    if (!ibCellsPtr_)
    {
        calcImmersedBoundary();
    }

    return *ibCellsPtr_;
}


const Foam::vectorField&
Foam::immersedBoundaryPolyPatch::ibCellCentres() const
{
    if (!ibCellCentresPtr_)
    {
        calcImmersedBoundary();
    }

    return *ibCellCentresPtr_;
}


const Foam::scalarField&
Foam::immersedBoundaryPolyPatch::ibCellVolumes() const
{
    if (!ibCellVolumesPtr_)
    {
        calcImmersedBoundary();
    }

    return *ibCellVolumesPtr_;
}


const Foam::labelList&
Foam::immersedBoundaryPolyPatch::ibFaces() const
{
    if (!ibFacesPtr_)
    {
        calcImmersedBoundary();
    }

    return *ibFacesPtr_;
}


const Foam::vectorField&
Foam::immersedBoundaryPolyPatch::ibFaceCentres() const
{
    if (!ibFaceCentresPtr_)
    {
        calcImmersedBoundary();
    }

    return *ibFaceCentresPtr_;
}


const Foam::vectorField&
Foam::immersedBoundaryPolyPatch::ibFaceAreas() const
{
    if (!ibFaceAreasPtr_)
    {
        calcImmersedBoundary();
    }

    return *ibFaceAreasPtr_;
}


const Foam::labelList&
Foam::immersedBoundaryPolyPatch::nearestTri() const
{
    if (!nearestTriPtr_)
    {
        calcImmersedBoundary();
    }

    return *nearestTriPtr_;
}


const Foam::labelList&
Foam::immersedBoundaryPolyPatch::deadCells() const
{
    if (!deadCellsPtr_)
    {
        calcImmersedBoundary();
    }

    return *deadCellsPtr_;
}


const Foam::labelList&
Foam::immersedBoundaryPolyPatch::deadFaces() const
{
    if (!deadFacesPtr_)
    {
        calcImmersedBoundary();
    }

    return *deadFacesPtr_;
}


const Foam::vectorField&
Foam::immersedBoundaryPolyPatch::correctedCellCentres() const
{
    if (!correctedIbPatchFaceAreasPtr_)
    {
        calcCorrectedGeometry();
    }

    return correctedFields_.correctedCellCentres();
}


const Foam::vectorField&
Foam::immersedBoundaryPolyPatch::correctedFaceCentres() const
{
    if (!correctedIbPatchFaceAreasPtr_)
    {
        calcCorrectedGeometry();
    }

    return correctedFields_.correctedFaceCentres();
}


const Foam::scalarField&
Foam::immersedBoundaryPolyPatch::correctedCellVolumes() const
{
    if (!correctedIbPatchFaceAreasPtr_)
    {
        calcCorrectedGeometry();
    }

    return correctedFields_.correctedCellVolumes();
}


const Foam::vectorField&
Foam::immersedBoundaryPolyPatch::correctedFaceAreas() const
{
    if (!correctedIbPatchFaceAreasPtr_)
    {
        calcCorrectedGeometry();
    }

    return correctedFields_.correctedFaceAreas();
}


const Foam::vectorField&
Foam::immersedBoundaryPolyPatch::correctedIbPatchFaceAreas() const
{
    if (!correctedIbPatchFaceAreasPtr_)
    {
        calcCorrectedGeometry();
    }

    return *correctedIbPatchFaceAreasPtr_;
}


const Foam::pointField&
Foam::immersedBoundaryPolyPatch::oldIbPoints() const
{
    if (!oldIbPointsPtr_)
    {
        // The mesh has never moved: old points are equal to current points
        ibUpdateTimeIndex_ = boundaryMesh().mesh().time().timeIndex();

        oldIbPointsPtr_ = new pointField(ibMesh_.points());
    }

    return *oldIbPointsPtr_;
}

Foam::tmp<Foam::vectorField>
Foam::immersedBoundaryPolyPatch::triMotionDistance() const
{
    // Calculate the distance between new and old coordinates on
    // the ibPatch face centres

    // Calculate the motion on the triangular mesh face centres
    return ibMesh_.coordinates()
      - PrimitivePatch<labelledTri, List, const pointField&>
        (
            ibMesh_,
            oldIbPoints()
        ).faceCentres();
}


Foam::tmp<Foam::vectorField>
Foam::immersedBoundaryPolyPatch::motionDistance() const
{
    // Interpolate the values from tri surface using nearest triangle
    return tmp<vectorField>
    (
        new vectorField(triMotionDistance(), nearestTri())
    );
}


void Foam::immersedBoundaryPolyPatch::moveTriSurfacePoints
(
    const pointField& p
)
{
    // Record the motion of the patch
    movingIb_ = true;

    // Move points of the triSurface
    const pointField& oldPoints = ibMesh_.points();

    if (oldPoints.size() != p.size())
    {
        FatalErrorIn
        (
            "void immersedBoundaryPolyPatch::moveTriSurfacePoints\n"
            "(\n"
            "    const pointField& p\n"
            ")"
        )   << "Incorrect size of motion points for patch " << name()
            << ".  oldPoints = "
            << oldPoints.size() << " p = " << p.size()
            << abort(FatalError);
    }

    if (ibUpdateTimeIndex_ < boundaryMesh().mesh().time().timeIndex())
    {
        // New motion in the current time step.  Store old points
        ibUpdateTimeIndex_ = boundaryMesh().mesh().time().timeIndex();

        deleteDemandDrivenData(oldIbPointsPtr_);
        Info<< "Storing old points for time index " << ibUpdateTimeIndex_
            << endl;
        oldIbPointsPtr_ = new pointField(oldPoints);
    }

    Info<< "Moving immersed boundary points for patch " << name()
        << endl;

    ibMesh_.movePoints(p);

    if (boundaryMesh().mesh().time().outputTime())
    {
        fileName path(boundaryMesh().mesh().time().path()/"VTK");

        mkDir(path);
        ibMesh_.triSurface::write
        (
            path/
            word
            (
                name() + "_tri_"
                + Foam::name(boundaryMesh().mesh().time().timeIndex())
                + ".stl"
            )
        );
    }

    // Note: the IB patch is now in the new position, but the mesh has not
    // been updated yet.  movePoints() needs to be executed to update the
    // fv mesh data
}


void Foam::immersedBoundaryPolyPatch::clearGeom()
{
    clearOut();
}


void Foam::immersedBoundaryPolyPatch::clearAddressing()
{
    clearOut();
}


void Foam::immersedBoundaryPolyPatch::clearOut() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryPolyPatch::clearOut() const")
            << "Clear immersed boundary polyPatch"
            << endl;
    }

    deleteDemandDrivenData(triSurfSearchPtr_);

    deleteDemandDrivenData(ibPatchPtr_);
    deleteDemandDrivenData(ibCellsPtr_);
    deleteDemandDrivenData(ibCellCentresPtr_);
    deleteDemandDrivenData(ibCellVolumesPtr_);
    deleteDemandDrivenData(ibFacesPtr_);
    deleteDemandDrivenData(ibFaceCentresPtr_);
    deleteDemandDrivenData(ibFaceAreasPtr_);
    deleteDemandDrivenData(nearestTriPtr_);
    deleteDemandDrivenData(deadCellsPtr_);
    deleteDemandDrivenData(deadFacesPtr_);

    deleteDemandDrivenData(correctedIbPatchFaceAreasPtr_);

    // Update topo change index
    topoChangeIndex_++;

    // Clear global data MeshObject
    correctedFields_.clearOut(topoChangeIndex_);

    // Cannot delete old motion points.  HJ, 10/Dec/2017
}


void Foam::immersedBoundaryPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("internalFlow") << internalFlow_
        << token::END_STATEMENT << nl;
    os.writeKeyword("isWall") << isWall_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
