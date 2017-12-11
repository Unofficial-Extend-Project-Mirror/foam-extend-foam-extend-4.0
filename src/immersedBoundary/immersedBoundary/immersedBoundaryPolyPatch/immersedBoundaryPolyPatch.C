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
Foam::immersedBoundaryPolyPatch::liveFactor_
(
    "immersedBoundaryLiveFactor",
    1e-7
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::immersedBoundaryPolyPatch::cellSpan
(
    const label cellID
) const
{
    const polyMesh& mesh = boundaryMesh().mesh();

    // Calculate span as twice the bounding box size
    const scalar delta = 2*cmptMax
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
            // Intersected cell
            intersectedCell[cellI] = immersedPoly::CUT;
            nIntersectedCells++;
        }
    }

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
                span,
                internalFlow(),
                false               // iterate intersection
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
                // Real intesection.  Check cut
                const scalar cellFactor = cutCell.wetVolume()/V[cellI];

                if (cellFactor < liveFactor_())
                {
                    // Thin cut: cell is dry
                    Info<< "Dry cell from intersection" << cellI << endl;
                    intersectedCell[cellI] = immersedPoly::DRY;
                }
                else if (cellFactor > (1 - liveFactor_()))
                {
                    // Thick cut: cell is wet
                    Info<< "Wet cell from intersection" << cellI << endl;
                    intersectedCell[cellI] = immersedPoly::WET;
                }
                else
                {
                    // True intersection.  Cut the cell and store all
                    // derived data

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


    // Reset the cell lists
    Info<< "nIbCells: " << nIbCells << endl;
    unmergedFaces.setSize(nIbCells);
    ibCells.setSize(nIbCells);
    ibCellCentres.setSize(nIbCells);
    ibCellVolumes.setSize(nIbCells);
    nearestTri.setSize(nIbCells);

    // Check tri addressing
    if (min(nearestTri) == -1)
    {
        FatalErrorIn("immersedBoundaryPolyPatch::calcImmersedBoundary() const")
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

        Info << "Writing immersed patch as VTK" << endl;

        fileName fvPath(mesh.time().path()/"VTK");
        mkDir(fvPath);

        OStringStream outputFilename;
        outputFilename << "immersed" << name();

        ibPatchPtr_->writeVTK(fvPath/fileName(outputFilename.str()));

        outputFilename << "normals";
        ibPatchPtr_->writeVTKNormals(fvPath/fileName(outputFilename.str()));
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
        Info<< "nDeadCells: " << nDeadCells << endl;
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
        if
        (
            intersectedCell[owner[faceI]] == immersedPoly::WET
         && intersectedCell[neighbour[faceI]] == immersedPoly::WET
        )
        {
            intersectedFace[faceI] = immersedPoly::WET;
        }

        if
        (
            intersectedCell[owner[faceI]] == immersedPoly::DRY
         && intersectedCell[neighbour[faceI]] == immersedPoly::DRY
        )
        {
            intersectedFace[faceI] = immersedPoly::DRY;
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
    for (label faceI = mesh.nInternalFaces(); faceI < owner.size(); faceI++)
    {
        if
        (
            intersectedCell[owner[faceI]] == immersedPoly::WET
        )
        {
            intersectedFace[faceI] = immersedPoly::WET;
        }

        if
        (
            intersectedCell[owner[faceI]] == immersedPoly::DRY
        )
        {
            intersectedFace[faceI] = immersedPoly::DRY;
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
                    false               // iterate intersection
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
                    // Real intesection.  Check cut
                    const scalar faceFactor =
                        cutFace.wetAreaMag()/mag(S[faceI]);

                    if (faceFactor < liveFactor_())
                    {
                        // Thin cut: face is dry
                        intersectedFace[faceI] = immersedPoly::DRY;
                    }
                    else if (faceFactor > (1 - liveFactor_()))
                    {
                        // Thick cut: face is wet
                        intersectedFace[faceI] = immersedPoly::DRY;
                    }
                    else
                    {
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
    }

    // Reset the sizes of the list
    Info<< "nIbFaces: " << nIbFaces << endl;
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
        Info<< "nDeadFaces: " << nDeadFaces << endl;
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
    }
}


void Foam::immersedBoundaryPolyPatch::clearOut()
{
    deleteDemandDrivenData(triSurfSearchPtr_);
    deleteDemandDrivenData(ibPatchPtr_);
    deleteDemandDrivenData(ibCellsPtr_);
    deleteDemandDrivenData(ibCellCentresPtr_);
    deleteDemandDrivenData(ibFacesPtr_);
    deleteDemandDrivenData(ibFaceCentresPtr_);
    deleteDemandDrivenData(nearestTriPtr_);
    deleteDemandDrivenData(deadCellsPtr_);
    deleteDemandDrivenData(deadFacesPtr_);

    deleteDemandDrivenData(correctedCellCentresPtr_);
    deleteDemandDrivenData(correctedFaceCentresPtr_);
    deleteDemandDrivenData(correctedCellVolumesPtr_);
    deleteDemandDrivenData(correctedFaceAreasPtr_);
    deleteDemandDrivenData(correctedIbPatchFaceAreasPtr_);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

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
        correctedCellCentresPtr_
     || correctedFaceCentresPtr_
     || correctedCellVolumesPtr_
     || correctedFaceAreasPtr_
     || correctedIbPatchFaceAreasPtr_
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

    // Copy basic fields for corrections
    correctedCellCentresPtr_ = new vectorField(mesh.cellCentres());
    vectorField& C = *correctedCellCentresPtr_;

    correctedFaceCentresPtr_ = new vectorField(mesh.faceCentres());
    vectorField& Cf = *correctedFaceCentresPtr_;

    correctedCellVolumesPtr_ = new scalarField(mesh.cellVolumes());
    scalarField& V = *correctedCellVolumesPtr_;

    correctedFaceAreasPtr_ = new vectorField(mesh.faceAreas());
    vectorField& Sf = *correctedFaceAreasPtr_;

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
        // Set dead volume to small
        V[dc[dcI]] = SMALL;
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
        // Set dead area to small
        Sf[df[dfI]] *= SMALL;
    }

    // In case of cutting errors due to finite tolerance, some cut cells may
    // remain opened and have to be closed by force.  This will be achieved
    // by the Marooney Maneouvre, where the face sum imbalance is compensated
    // in the cut face.  HJ, 11/Dec/2017

    vectorField sumSf(mesh.nCells(), vector::zero);

    // Get face addressing
    const labelList& owner = mesh.faceOwner();
    const labelList& neighbour = mesh.faceNeighbour();

    forAll (owner, faceI)
    {
        sumSf[owner[faceI]] += Sf[faceI];
    }

    forAll (neighbour, faceI)
    {
        sumSf[neighbour[faceI]] -= Sf[faceI];
    }

    forAll (cutCells, cutCellI)
    {
        const label ccc = cutCells[cutCellI];
        if
        (
            mag(sumSf[ccc]) > 1e-12
         && mag(sumSf[ccc] + ibSf[cutCellI])/cutCellVolumes[cutCellI] > 1e-6
        )
        {
            Info<< "Marooney Maneouvre for cell " << ccc
                << " error: " << mag(sumSf[ccc] + ibSf[cutCellI]) << " "
                << mag(sumSf[ccc] + ibSf[cutCellI])/cutCellVolumes[cutCellI]
                << " " << sumSf[ccc] << " "
                << " V: " << cutCellVolumes[cutCellI]
                << " Sf: " << ibSf[cutCellI] << endl;
            ibSf[cutCellI] = -sumSf[ccc];
        }
    }
}


void Foam::immersedBoundaryPolyPatch::movePoints(const pointField& p)
{
    // Handle motion of immersed boundary
    clearOut();

    primitivePatch::movePoints(p);
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
            bm.mesh(),                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    internalFlow_(false),
    isWall_(true),
    movingIb_(false),
    ibUpdateTimeIndex_(-1),
    triSurfSearchPtr_(NULL),
    ibPatchPtr_(NULL),
    ibCellsPtr_(NULL),
    ibCellCentresPtr_(NULL),
    ibFacesPtr_(NULL),
    ibFaceCentresPtr_(NULL),
    nearestTriPtr_(NULL),
    deadCellsPtr_(NULL),
    deadFacesPtr_(NULL),
    correctedCellCentresPtr_(NULL),
    correctedFaceCentresPtr_(NULL),
    correctedCellVolumesPtr_(NULL),
    correctedFaceAreasPtr_(NULL),
    correctedIbPatchFaceAreasPtr_(NULL)
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
            bm.mesh(),                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    internalFlow_(dict.lookup("internalFlow")),
    isWall_(dict.lookup("isWall")),
    movingIb_(false),
    ibUpdateTimeIndex_(-1),
    triSurfSearchPtr_(NULL),
    ibPatchPtr_(NULL),
    ibCellsPtr_(NULL),
    ibCellCentresPtr_(NULL),
    ibCellVolumesPtr_(NULL),
    ibFacesPtr_(NULL),
    ibFaceCentresPtr_(NULL),
    ibFaceAreasPtr_(NULL),
    nearestTriPtr_(NULL),
    deadCellsPtr_(NULL),
    deadFacesPtr_(NULL),
    correctedCellCentresPtr_(NULL),
    correctedFaceCentresPtr_(NULL),
    correctedCellVolumesPtr_(NULL),
    correctedFaceAreasPtr_(NULL),
    correctedIbPatchFaceAreasPtr_(NULL)
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
            bm.mesh(),                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    internalFlow_(pp.internalFlow_),
    isWall_(pp.isWall_),
    movingIb_(false),
    ibUpdateTimeIndex_(-1),
    triSurfSearchPtr_(NULL),
    ibPatchPtr_(NULL),
    ibCellsPtr_(NULL),
    ibCellCentresPtr_(NULL),
    ibFacesPtr_(NULL),
    ibFaceCentresPtr_(NULL),
    nearestTriPtr_(NULL),
    deadCellsPtr_(NULL),
    deadFacesPtr_(NULL),
    correctedCellCentresPtr_(NULL),
    correctedFaceCentresPtr_(NULL),
    correctedCellVolumesPtr_(NULL),
    correctedFaceAreasPtr_(NULL),
    correctedIbPatchFaceAreasPtr_(NULL)
{}


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
            bm.mesh(),                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    internalFlow_(pp.internalFlow_),
    isWall_(pp.isWall_),
    movingIb_(false),
    ibUpdateTimeIndex_(-1),
    triSurfSearchPtr_(NULL),
    ibPatchPtr_(NULL),
    ibCellsPtr_(NULL),
    ibCellCentresPtr_(NULL),
    ibFacesPtr_(NULL),
    ibFaceCentresPtr_(NULL),
    nearestTriPtr_(NULL),
    deadCellsPtr_(NULL),
    deadFacesPtr_(NULL),
    correctedCellCentresPtr_(NULL),
    correctedFaceCentresPtr_(NULL),
    correctedCellVolumesPtr_(NULL),
    correctedFaceAreasPtr_(NULL),
    correctedIbPatchFaceAreasPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::immersedBoundaryPolyPatch::~immersedBoundaryPolyPatch()
{
    clearOut();
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
    if (!correctedCellCentresPtr_)
    {
        calcCorrectedGeometry();
    }

    return *correctedCellCentresPtr_;
}


const Foam::vectorField&
Foam::immersedBoundaryPolyPatch::correctedFaceCentres() const
{
    if (!correctedFaceCentresPtr_)
    {
        calcCorrectedGeometry();
    }

    return *correctedFaceCentresPtr_;
}


const Foam::scalarField&
Foam::immersedBoundaryPolyPatch::correctedCellVolumes() const
{
    if (!correctedCellVolumesPtr_)
    {
        calcCorrectedGeometry();
    }

    return *correctedCellVolumesPtr_;
}


const Foam::vectorField&
Foam::immersedBoundaryPolyPatch::correctedFaceAreas() const
{
    if (!correctedFaceAreasPtr_)
    {
        calcCorrectedGeometry();
    }

    return *correctedFaceAreasPtr_;
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

    Info<< "Moving immersed boundary points for patch " << name()
        << endl;

    ibMesh_.movePoints(p);

    fileName path(boundaryMesh().mesh().time().path()/"VTK");

    mkDir(path);
    ibMesh_.triSurface::write
    (
        path/
        word
        (
            name() + "_"
          + Foam::name(boundaryMesh().mesh().time().timeIndex())
          + ".stl"
        )
    );

    clearOut();
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
