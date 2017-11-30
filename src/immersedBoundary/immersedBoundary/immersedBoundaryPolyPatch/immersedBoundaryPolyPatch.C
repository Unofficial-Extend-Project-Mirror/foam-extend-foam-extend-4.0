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
    1e-9
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::immersedBoundaryPolyPatch::cellSpan
(
    const label cellID
) const
{
    const polyMesh& mesh = boundaryMesh().mesh();

    scalar delta;

    if (mesh.nGeometricD() == 3)
    {
        delta = 2*pow(mesh.cellVolumes()[cellID], 1.0/3.0);
    }
    else
    {
        const Vector<label>& directions = mesh.geometricD();

        scalar thickness = 0.0;

        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh.bounds().span()[dir];
                break;
            }
        }

        // Field created with mapping for IB cells only
        delta = 2*sqrt(mesh.cellVolumes()[cellID]/thickness);
    }

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
     || gammaPtr_
     || sGammaPtr_
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

    // Get mesh face areas
    const vectorField& S = mesh.faceAreas();

    // Get mesh cell volumes
    const scalarField& V = mesh.cellVolumes();

    // Get mesh cell centres
    const vectorField& C = mesh.cellCentres();

    // Get face addressing
    const labelList& owner = mesh.faceOwner();
    const labelList& neighbour = mesh.faceNeighbour();

    // Algorithm
    // Initialise the search by marking the inside cells using calcInside
    // Based on inside cells addressing, for every mesh face straddling
    // the surface, check if the cell straddles the free surface
    // For cells next to the coupled boundary, check if the points on the
    // boundary have a different inside index from the cell centres next to it.
    // Collect all intersected cells
    // For all intersected cells, calculate the actual intersection and
    // - calculate the intersection face, its centre, and area vector
    // - adjust the cell volume, using a cell blending factor 0 < gamma < 1
    // - adjust the face area for all cell faces, using a
    //   face blending factor 0 < gamma < 1

    // Mark cells that are inside or outside of the triangular surface
    boolList cellsInside = tss.calcInside(C);

    // Initialise intersected
    boolList intersectedCell(mesh.nCells(), false);

    // Go through the faces at the interface between a live and dead cell
    // and mark the band of possible intersections
    forAll (neighbour, faceI)
    {
        if (cellsInside[owner[faceI]] != cellsInside[neighbour[faceI]])
        {
            // Check both cell for intersection
            intersectedCell[owner[faceI]] = true;
            intersectedCell[neighbour[faceI]] = true;
        }
    }

    // Go through the patches and repeat the check with patch points
    // Note: checking all patches for IB cut close to the boundary
    // HJ. 24/Nov/2017
    forAll (bMesh, patchI)
    {
        // Check all cells next to a patch
        const polyPatch& curPatch = bMesh[patchI];

        // Skip empty patch
        if (isA<emptyPolyPatch>(curPatch))
        {
            continue;
        }

        const labelList& faceCells = curPatch.faceCells();
        const faceList& localFaces = curPatch.localFaces();

        // For all points on the patch, check intersection status
        boolList facePointInside = tss.calcInside(curPatch.localPoints());

        // For each cell next to the patch, check whether all points
        // are on the same side of the patch as the cell centre

        forAll (faceCells, faceI)
        {
            // Get cell centre status
            const bool curCellInside = cellsInside[faceCells[faceI]];

            const face& curFace = localFaces[faceI];

            forAll (curFace, pointI)
            {
                if (facePointInside[curFace[pointI]] != curCellInside)
                {
                    // Cell is cut.  Mark it
                    intersectedCell[faceCells[faceI]] = true;

                    break;
                }
            }
        }
    }

    // Count intersected cells and allocate memory
    // Note: some cells marked as possibly intersected may be rejected later
    label nIbCells = 0;

    forAll (intersectedCell, cellI)
    {
        if (intersectedCell[cellI])
        {
            nIbCells++;
        }
    }
    Info<< "nIntersectedCells: " << nIbCells << endl;
    // Collect intersection points and faces.  Primitive patch will be created
    // after renumbering

    // IB points
    DynamicList<point> unmergedPoints(nIbCells*primitiveMesh::pointsPerFace_);
    label nIbPoints = 0;

    // IB patch faces: Cell intersections with the IB patch
    faceList unmergedFaces(nIbCells);

    // IB cells: cells intersected by the IB patch
    // This also corresponds to faceCells next to the IB patch
    ibCellsPtr_ = new labelList(nIbCells);
    labelList& ibCells = *ibCellsPtr_;

    // IB cellCentres: centre of live part of the intersected cell
    // next to the IB patch
    ibCellCentresPtr_ = new vectorField(nIbCells);
    vectorField& ibCellCentres = *ibCellCentresPtr_;

    // IB cellCentres: centre of live part of the intersected cell
    // next to the IB patch
    ibCellVolumesPtr_ = new scalarField(nIbCells);
    scalarField& ibCellVolumes = *ibCellVolumesPtr_;

    // Nearest triangle
    nearestTriPtr_ = new labelList(nIbCells);
    labelList& nearestTri = *nearestTriPtr_;

    // Reset the counter
    nIbCells = 0;

    // Calculate gamma
    gammaPtr_ = new scalarField(mesh.nCells(), 1);
    scalarField& gamma = *gammaPtr_;

    // Collect dead cells
    boolList deadCells(mesh.nCells(), false);

    // Adjust selection of cells: inside or outside of immersed boundary
    if (internalFlow())
    {
        Info<< "Internal flow" << endl;
        forAll (gamma, cellI)
        {
            if (!cellsInside[cellI])
            {
                gamma[cellI] = 0;
                deadCells[cellI] = true;
            }
        }
    }
    else
    {
        Info<< "External flow" << endl;
        forAll (gamma, cellI)
        {
            if (cellsInside[cellI])
            {
                gamma[cellI] = 0;
                deadCells[cellI] = true;
            }
        }
    }

    // At this point, all live cells are marked with 1
    // Intesect all cells that are marked for intersection

    forAll (intersectedCell, cellI)
    {
        if (intersectedCell[cellI])
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
                gamma[cellI] = 1;
            }
            else if (cutCell.isAllDry())
            {
                gamma[cellI] = 0;
                deadCells[cellI] = true;

            }
            else
            {
                // Real intesection.  Check cut
                const scalar cellFactor = cutCell.wetVolume()/V[cellI];

                if (cellFactor < liveFactor_())
                {
                    // Thin cut: cell is dry
                    gamma[cellI] = 0;
                    deadCells[cellI] = true;
                }
                else if (cellFactor > (1 - liveFactor_()))
                {
                    // Thick cut: cell is wet
                    gamma[cellI] = 0;
                    deadCells[cellI] = true;
                }
                else
                {
                    // True intersection.  Cut the cell and store all
                    // derived data
                    deadCells[cellI] = false;

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

                    // Record the live volume
                    gamma[cellI] = cellFactor;
                }
            }
        }
    }

    // Reset the cell lists
    Info<< "nIbCells: " << nIbCells << endl;
    unmergedFaces.setSize(nIbCells);
    ibCells.setSize(nIbCells);
    ibCellCentres.setSize(nIbCells);
    ibCellVolumes.setSize(nIbCells);

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

        forAll (deadCells, cellI)
        {
            if (deadCells[cellI])
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

        forAll (deadCells, cellI)
        {
            if (deadCells[cellI])
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

    // Calculate sGamma
    sGammaPtr_ = new scalarField(mesh.nFaces(), 1);
    scalarField& sGamma = *sGammaPtr_;

    // Collect dead faces
    boolList deadFaces(mesh.nFaces(), false);

    // First, kill all faces touching dead cells, including internal
    // and boundary faces
    // The intersection belt will be handled separately
    forAll (neighbour, faceI)
    {
        if (deadCells[neighbour[faceI]])
        {
            deadFaces[faceI] = true;
        }
    }

    // Internal and boundary faces
    forAll (owner, faceI)
    {
        if (deadCells[owner[faceI]])
        {
            deadFaces[faceI] = true;
        }
    }
    
    // For all internal faces, check if the owner or neighbour cell has been cut
    // For all boundary faces, check if the internal cell has been cut

    // Internal face check
    forAll (neighbour, faceI)
    {
        if
        (
            intersectedCell[owner[faceI]]
         || intersectedCell[neighbour[faceI]]
        )
        {
            // Possibly intersected face  Try to cut it
            // Create a cutting object with a local tolerance
            triSurfaceDistance dist
            (
                tss,
                Foam::max
                (
                    cellSpan(owner[faceI]),
                    cellSpan(neighbour[faceI])
                ),
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
                sGamma[faceI] = 1;
            }
            else if (cutFace.isAllDry())
            {
                sGamma[faceI] = 0;
                deadFaces[faceI] = true;
            }
            else
            {
                // Real intesection.  Check cut
                const scalar faceFactor = cutFace.wetAreaMag()/mag(S[faceI]);

                if (faceFactor < liveFactor_())
                {
                    // Thin cut: face is dry
                    sGamma[faceI] = 0;
                    deadFaces[faceI] = true;
                }
                else if (faceFactor > (1 - liveFactor_()))
                {
                    // Thick cut: face is wet
                    sGamma[faceI] = 0;
                }
                else
                {
                    // True intersection.  Collect data

                    // Get intersected face index
                    ibFaces[nIbFaces] = faceI;

                    // Get wet centre
                    ibFaceCentres[nIbFaces] = cutFace.wetAreaCentre();

                    // Get wet area
                    ibFaceAreas[nIbFaces] = faceFactor*S[faceI];

                    // Calculate wet fraction from wet area magnitude and
                    // original area magnitude
                    sGamma[faceI] = faceFactor;

                    nIbFaces++;
                }
            }
        }
        else
        {
            // No intersection

            // Debug check: if neither owner nor neighbour has been cut, they
            // need to have the same value of gamma
            // Comparing 0:1 indicator function
            if (mag(gamma[owner[faceI]] - gamma[neighbour[faceI]]) > SMALL)
            {
                FatalErrorIn
                (
                    "void immersedBoundaryPolyPatch::"
                    "calcImmersedBoundary() const"
                )
                    << "Topological face cutting error for patch "
                    << name() << ", face " << faceI
                    << abort(FatalError);
            }

            sGamma[faceI] = gamma[owner[faceI]];
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

        forAll (deadFaces, faceI)
        {
            if (deadFaces[faceI])
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

        forAll (deadFaces, faceI)
        {
            if (deadFaces[faceI])
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
    deleteDemandDrivenData(gammaPtr_);
    deleteDemandDrivenData(sGammaPtr_);

    deleteDemandDrivenData(correctedCellCentresPtr_);
    deleteDemandDrivenData(correctedFaceCentresPtr_);
    deleteDemandDrivenData(correctedCellVolumesPtr_);
    deleteDemandDrivenData(correctedFaceAreasPtr_);
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

    // Correct for all cut cells

    // Memory management
    {
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
    }

    // Correct for all cut faces
    // Memory management
    {
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
    gammaPtr_(NULL),
    sGammaPtr_(NULL),
    correctedCellCentresPtr_(NULL),
    correctedFaceCentresPtr_(NULL),
    correctedCellVolumesPtr_(NULL),
    correctedFaceAreasPtr_(NULL)
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
    gammaPtr_(NULL),
    sGammaPtr_(NULL),
    correctedCellCentresPtr_(NULL),
    correctedFaceCentresPtr_(NULL),
    correctedCellVolumesPtr_(NULL),
    correctedFaceAreasPtr_(NULL)
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
    gammaPtr_(NULL),
    sGammaPtr_(NULL),
    correctedCellCentresPtr_(NULL),
    correctedFaceCentresPtr_(NULL),
    correctedCellVolumesPtr_(NULL),
    correctedFaceAreasPtr_(NULL)
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
    gammaPtr_(NULL),
    sGammaPtr_(NULL),
    correctedCellCentresPtr_(NULL),
    correctedFaceCentresPtr_(NULL),
    correctedCellVolumesPtr_(NULL),
    correctedFaceAreasPtr_(NULL)
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


const Foam::scalarField&
Foam::immersedBoundaryPolyPatch::gamma() const
{
    if (!gammaPtr_)
    {
        calcImmersedBoundary();
    }

    return *gammaPtr_;
}


const Foam::scalarField&
Foam::immersedBoundaryPolyPatch::sGamma() const
{
    if (!sGammaPtr_)
    {
        calcImmersedBoundary();
    }

    return *sGammaPtr_;
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
}


// ************************************************************************* //
