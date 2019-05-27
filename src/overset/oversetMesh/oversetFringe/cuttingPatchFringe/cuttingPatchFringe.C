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

#include "cuttingPatchFringe.H"
#include "oversetMesh.H"
#include "oversetRegion.H"
#include "polyPatchID.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cuttingPatchFringe, 0);
    addToRunTimeSelectionTable
    (
        oversetFringe,
        cuttingPatchFringe,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cuttingPatchFringe::calcAddressing() const
{
    // Make sure that either acceptorsPtr is unnalocated or if it is allocated,
    // that it is empty
    if (acceptorsPtr_ && !acceptorsPtr_->empty())
    {
        FatalErrorIn
        (
            "void Foam::cuttingPatchFringe::calcAddressing() const"
        )   << "Addressing already calculated"
            << abort(FatalError);
    }

    // Get polyMesh
    const polyMesh& mesh = this->mesh();

    // Collect all cutting patches
    labelHashSet patchIDs(cuttingPatchNames_.size());

    forAll (cuttingPatchNames_, nameI)
    {
        // Get polyPatchID and check if valid
        const polyPatchID cutPatch
        (
            cuttingPatchNames_[nameI],
            mesh.boundaryMesh()
        );

        if (!cutPatch.active())
        {
            FatalErrorIn
            (
                "void cuttingPatchFringe::calcAddressing const"
            )   << "Cutting patch " << cuttingPatchNames_[nameI]
                << " cannot be found."
                << abort(FatalError);
        }

        // Store patch ID in the set
        patchIDs.insert(cutPatch.index());
    }

    if (debug)
    {
        Info<< "Starting cutting patch fringe assembly..." << endl;
    }

    // Note: similar code as in oversetRegion::calcHoleTriMesh. Consider
    // refactoring. VV, 20/May/2019

    // Make and invert local triSurface
    triFaceList triFaces;
    pointField triPoints;

    // Memory management
    {
        triSurface ts = triSurfaceTools::triangulate
        (
            mesh.boundaryMesh(),
            patchIDs
        );

        // Clean mutiple points and zero-sized triangles
        ts.cleanup(false);

        triFaces.setSize(ts.size());
        triPoints = ts.points();

        forAll (ts, tsI)
        {
            // Bugfix: no need to reverse the face because the normals point in
            // the correct direction already. VV, 20/May/2019.
            triFaces[tsI] = ts[tsI];
        }
    }

    if (Pstream::parRun())
    {
        // Combine all faces and points into a single list

        List<triFaceList> allTriFaces(Pstream::nProcs());
        List<pointField> allTriPoints(Pstream::nProcs());

        allTriFaces[Pstream::myProcNo()] = triFaces;
        allTriPoints[Pstream::myProcNo()] = triPoints;

        Pstream::gatherList(allTriFaces);
        Pstream::scatterList(allTriFaces);

        Pstream::gatherList(allTriPoints);
        Pstream::scatterList(allTriPoints);

        // Re-pack points and faces

        label nTris = 0;
        label nPoints = 0;

        forAll (allTriFaces, procI)
        {
            nTris += allTriFaces[procI].size();
            nPoints += allTriPoints[procI].size();
        }

        // Pack points
        triPoints.setSize(nPoints);

        // Prepare point renumbering array
        labelListList renumberPoints(Pstream::nProcs());

        nPoints = 0;

        forAll (allTriPoints, procI)
        {
            const pointField& ptp = allTriPoints[procI];

            renumberPoints[procI].setSize(ptp.size());

            labelList& procRenumberPoints = renumberPoints[procI];

            forAll (ptp, ptpI)
            {
                triPoints[nPoints] = ptp[ptpI];
                procRenumberPoints[ptpI] = nPoints;

                nPoints++;
            }
        }

        // Pack triangles and renumber into complete points on the fly
        triFaces.setSize(nTris);

        nTris = 0;

        forAll (allTriFaces, procI)
        {
            const triFaceList& ptf = allTriFaces[procI];

            const labelList& procRenumberPoints = renumberPoints[procI];

            forAll (ptf, ptfI)
            {
                const triFace& procFace = ptf[ptfI];

                triFace& renumberFace = triFaces[nTris];

                forAll (renumberFace, rfI)
                {
                    renumberFace[rfI] =
                        procRenumberPoints[procFace[rfI]];
                }

                nTris++;
            }
        }
    }

    // Make a complete triSurface from local data
    triSurface patchTriMesh(triFaces, triPoints);

    // Clean up duplicate points and zero sized triangles
    patchTriMesh.cleanup(false);

    // Get this region
    const oversetRegion& myRegion = this->region();

    // Debug: write down the tri mesh
    if (Pstream::master())
    {
        patchTriMesh.write
        (
            word
            (
                "patchTriSurface_region" + myRegion.name() + ".vtk"
            )
        );
    }

    // Create the tri surface search object
    triSurfaceSearch patchTriSearch(patchTriMesh);

    // Get cells in this region
    const labelList& myRegionCells = myRegion.regionCells();

    // Get cell centres for inside-outside search using search object
    vectorField myCC(myRegionCells.size());

    // Cell centres from polyMesh
    const vectorField& cc = mesh.cellCentres();

    forAll (myCC, i)
    {
        myCC[i] = cc[myRegionCells[i]];
    }

    // Inside mask: all cells within search object will be marked
    boolList insideMask(mesh.nCells(), false);

    // Get inside cells for cells in my region only
    boolList myRegionInsideMask = patchTriSearch.calcInside(myCC);

    // Note: insideMask has the size of all mesh cells and
    // myRegionInsideMask has the size of cells in this region
    forAll (myRegionInsideMask, i)
    {
        insideMask[myRegionCells[i]] = myRegionInsideMask[i];
    }

    // Make sure that the cut holes for this region are also properly marked as
    // "inside". This may not be the case automatically for e.g. simulations
    // with appendages
    const labelList& cutRegionHoles = myRegion.cutHoles();
    forAll (cutRegionHoles, i)
    {
        insideMask[cutRegionHoles[i]] = true;
    }

    // Get necessary mesh data (from polyMesh/primitiveMesh)
    const cellList& meshCells = mesh.cells();
    const unallocLabelList& owner = mesh.faceOwner();
    const unallocLabelList& neighbour = mesh.faceNeighbour();

    // Bool list for collecting faces with at least one unmarked
    // cell (to determine the acceptors for the first iteration)
    boolList hasUnmarkedCell(mesh.nFaces(), false);

    // Loop through all cells
    forAll (insideMask, cellI)
    {
        if (!insideMask[cellI])
        {
            // This cell is not inside (it is unmarked). Loop through
            // its faces and set the flag
            const cell& cFaces = meshCells[cellI];

            forAll (cFaces, i)
            {
                // Set the mark for this global face
                hasUnmarkedCell[cFaces[i]] = true;
            }
        }
    }

    // Sync the face list across processor boundaries
    syncTools::syncFaceList
    (
        mesh,
        hasUnmarkedCell,
        orEqOp<bool>(),
        true
    );

    // Mark-up for all inside faces
    boolList insideFaceMask(mesh.nFaces(), false);

    // Collect all acceptors for the first iteration (the cells that
    // have at least one neighbour cell that is not marked)
    labelHashSet acceptors(myRegionCells.size()/10);

    // Loop again through all cells and collect marked ones into
    // acceptors or holes, depending on whether they have unmarked cell
    // as a neighbour (indicating an acceptor)
    forAll (insideMask, cellI)
    {
        if (insideMask[cellI])
        {
            // This cell is inside the covered region
            const cell& cFaces = meshCells[cellI];

            forAll (cFaces, i)
            {
                // Get global face index
                const label& faceI = cFaces[i];

                // Check whether this neighbour is unmarked
                if (hasUnmarkedCell[faceI])
                {
                    // This cell has unmarked neighbour, collect it into
                    // the acceptor list
                    acceptors.insert(cellI);

                    // This cell is no longer "inside cell"
                    insideMask[cellI] = false;
                }
                else
                {
                    // This is an "inside" face, mark it
                    insideFaceMask[faceI] = true;
                }
            } // End for all faces
        } // End if cell is inside
    } // End for all cells

    // Note: insideFaceMask already synced across processors because it relies
    // on hasUnmarkedCell list, which has been synced just above

    // Hash set containing new acceptors (for successive iterations)
    labelHashSet newAcceptors(acceptors.size());

    // Now that we have the initial set of acceptors (and holes), loop
    // nLayers away from initial donors
    for (label i = 0; i < nLayers_; ++i)
    {
        // Face markup for propagation
        boolList propagateFace(mesh.nFaces(), false);

        // Loop through all acceptors and mark faces that are around hole
        // cells. This way, we make sure that we go towards the correct,
        // inside direction
        forAllConstIter (labelHashSet, acceptors, iter)
        {
            // Get the cell index and the cell
            const label& cellI = iter.key();
            const cell& cFaces = meshCells[cellI];

            // Loop through all faces of the cell
            forAll (cFaces, i)
            {
                // Get face index (global)
                const label& faceI = cFaces[i];

                if (insideFaceMask[faceI])
                {
                    // This is a hole face, we are moving in the right
                    // direction. Mark the face for propagation
                    propagateFace[faceI] = true;
                }
            } // End for all faces of the cell
        } // End for all donor cells

        // Sync the face list across processor boundaries
        syncTools::syncFaceList
        (
            mesh,
            propagateFace,
            orOp<bool>(),
            false
        );

        // Loop through all internal faces and append acceptors
        for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
        {
            if (propagateFace[faceI])
            {
                // Face is marked, select owner or neighbour
                const label& own = owner[faceI];
                const label& nei = neighbour[faceI];

                // Either owner or neighbour may be hole, not both
                if (insideMask[own])
                {
                    // Owner cell is a hole, insert it
                    newAcceptors.insert(own);

                    // Update hole mask
                    insideMask[own] = false;
                }
                else if (insideMask[nei])
                {
                    // Neighbour cell is a hole, insert it
                    newAcceptors.insert(nei);

                    // Update hole mask
                    insideMask[nei] = false;
                }

                // Also update hole face mask for next iteration
                insideFaceMask[faceI] = false;
            }
        }

        // Loop through boundary faces
        for
        (
            label faceI = mesh.nInternalFaces();
            faceI < mesh.nFaces();
          ++faceI
        )
        {
            if (propagateFace[faceI])
            {
                // Face is marked, select owner if this is the right
                // side. Neighbour handled on the other side
                const label& own = owner[faceI];

                if (insideMask[own])
                {
                    // Face cell is a hole, insert it
                    newAcceptors.insert(own);

                    // Update hole mask
                    insideMask[own] = false;
                }

                // Also update hole face mask for next iteration
                insideFaceMask[faceI] = false;
            }
        }

        // Transfer newAcceptors into acceptors for next iteration or
        // for final assembly. Resize newAcceptors accordingly
        acceptors.transfer(newAcceptors);
        newAcceptors.resize(acceptors.size());

    } // End for specified number of layers

    // At this point, we have the final set of acceptors and we marked
    // all cells that should be holes. Collect them into the list
    dynamicLabelList fringeHoles(myRegionCells.size()/10);

    forAll (insideMask, cellI)
    {
        if (insideMask[cellI])
        {
            fringeHoles.append(cellI);
        }
    }

    // Set acceptors and holes from the data for all regions
    acceptorsPtr_ = new labelList(acceptors.sortedToc());
    fringeHolesPtr_ = new labelList(fringeHoles.xfer());

    if (debug)
    {
        Info<< "In cuttingPatchFringe::calcAddressing() const" << nl
            << "Found " << acceptorsPtr_->size() << " acceptors." << nl
            << "Found " << fringeHolesPtr_->size() << " fringe holes." << endl;
    }
}


void Foam::cuttingPatchFringe::clearAddressing() const
{
    deleteDemandDrivenData(fringeHolesPtr_);
    deleteDemandDrivenData(acceptorsPtr_);
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::cuttingPatchFringe::cuttingPatchFringe
(
    const fvMesh& mesh,
    const oversetRegion& region,
    const dictionary& dict
)
:
    oversetFringe(mesh, region, dict),
    cuttingPatchNames_(dict.lookup("cuttingPatches")),
    nLayers_(readLabel(dict.lookup("nLayers"))),
    fringeHolesPtr_(nullptr),
    acceptorsPtr_(nullptr),
    finalDonorAcceptorsPtr_(nullptr)
{
    // Sanity check number of layers: must be greater than 0
    if (nLayers_ < 1)
    {
        FatalIOErrorIn
        (
            "cuttingPatchFringe::"
            "cuttingPatchFringe\n"
            "(\n"
            "    const fvMesh& mesh,\n"
            "    const oversetRegion& region,\n"
            "    const dictionary& dict\n"
            ")",
            dict
        )   << "Invalid number of layers specified, nLayers = " << nLayers_
            << nl
            << "The number should be greater than 0."
            << abort(FatalError);
    }

    // Preferably, the number of layers should be at least 2
    if (nLayers_ == 1)
    {
        WarningIn
        (
            "cuttingPatchFringe::"
            "cuttingPatchFringe\n"
            "(\n"
            "    const fvMesh& mesh,\n"
            "    const oversetRegion& region,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "You have specified nLayers = " << nLayers_
            << nl
            << "We strongly advise to use at least 2 layers to avoid" << nl
            << "possibility of having acceptors that cannot find decent" << nl
            << "donors on the other side."
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cuttingPatchFringe::~cuttingPatchFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cuttingPatchFringe::updateIteration
(
    donorAcceptorList& donorAcceptorRegionData
) const
{
    // If the donorAcceptor list has been allocated, something went wrong with
    // the iteration procedure (not-updated flag): this function has been called
    // more than once, which should not happen for cuttingPatchFringe
    if (finalDonorAcceptorsPtr_)
    {
        FatalErrorIn
        (
            "cuttingPatchFringe::"
            "updateIteration(donorAcceptorList&) const"
        )   << "finalDonorAcceptorPtr_ already allocated. Something went "
            << "wrong with the iteration procedure (flag was not updated)."
            << nl
            << "This should not happen for cuttingPatchFringe."
            << abort(FatalError);
    }

    // Allocate the list by reusing the argument list
    finalDonorAcceptorsPtr_ = new donorAcceptorList
    (
        donorAcceptorRegionData,
        true
    );

    // Set the flag to true and return
    updateSuitableOverlapFlag(true);

    return foundSuitableOverlap();
}


const Foam::labelList& Foam::cuttingPatchFringe::fringeHoles() const
{
    if (!fringeHolesPtr_)
    {
        calcAddressing();
    }

    return *fringeHolesPtr_;
}


const Foam::labelList& Foam::cuttingPatchFringe::candidateAcceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }

    return *acceptorsPtr_;
}


Foam::donorAcceptorList&
Foam::cuttingPatchFringe::finalDonorAcceptors() const
{
    if (!finalDonorAcceptorsPtr_)
    {
        FatalErrorIn("cuttingPatchFringe::finalDonorAcceptors()")
            << "finalDonorAcceptorPtr_ not allocated. Make sure you have"
            << " called cuttingPatchFringe::updateIteration() before"
            << " asking for final set of donor/acceptor pairs."
            << abort(FatalError);
    }

    if (!foundSuitableOverlap())
    {
        FatalErrorIn("cuttingPatchFringe::finalDonorAcceptors()")
            << "Attemted to access finalDonorAcceptors but suitable overlap "
            << "has not been found. This is not allowed. "
            << abort(FatalError);
    }

    return *finalDonorAcceptorsPtr_;
}


void Foam::cuttingPatchFringe::update() const
{
    Info<< "cuttingPatchFringe::update() const" << endl;

    // Clear out
    clearAddressing();

    // Set flag to false and clear final donor/acceptors only
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
    updateSuitableOverlapFlag(false);
}


// ************************************************************************* //
