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

\*---------------------------------------------------------------------------*/

#include "oversetRegion.H"
#include "oversetMesh.H"
#include "oversetFringe.H"
#include "polyPatchID.H"
#include "triSurfaceTools.H"
#include "cellSet.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::oversetRegion::calcDonorRegions() const
{
    if (donorRegionsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcDonorRegions() const")
            << "Donor regions already calculated"
            << abort(FatalError);
    }

    // Get regions
    const PtrList<oversetRegion>& regions = oversetMesh_.regions();

    donorRegionsPtr_ = new labelList(donorRegionNames_.size());
    labelList& donRegions = *donorRegionsPtr_;

    forAll (donorRegionNames_, drI)
    {
        bool found = false;

        // Get name to search for
        const word& curName = donorRegionNames_[drI];

        // If the donor region name is the same as the name of this region,
        // issue an error
        if (name_ == curName)
        {
            FatalErrorIn
            (
                "void oversetRegion::calcDonorRegions() const"
            )   << "Region " << name_ << " specified as the donor "
                << "of itself.  List of donor regions: " << donorRegionNames_ << nl
                << "This is not allowed: check oversetMesh definition"
                << abort(FatalError);
        }

        // Find donor region name
        forAll (regions, orI)
        {
            if (regions[orI].name() == curName)
            {
                // Found donor region name in the list
                found = true;

                donRegions[drI] = orI;
                break;
            }
        }

        if (!found)
        {
            FatalErrorIn("void oversetRegion::calcDonorRegions() const")
                << "For region " << name() << " cannot find donor region "
                << curName << ".  Please check overset definition"
                << abort(FatalError);
        }
    }
}


void Foam::oversetRegion::calcAcceptorRegions() const
{
    if (acceptorRegionsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcAcceptorRegions() const")
            << "Acceptor regions already calculated"
            << abort(FatalError);
    }

    const PtrList<oversetRegion>& regions = oversetMesh_.regions();

    acceptorRegionsPtr_ = new labelList(regions.size());
    labelList& accRegions = *acceptorRegionsPtr_;

    label nAccRegions = 0;

    // Go through all regions apart from the current and check if
    // this region appears in the donor list
    forAll (regions, regionI)
    {
        // Skip current region
        if (regionI == index())
        {
            continue;
        }

        const labelList& remoteDonors = regions[regionI].donorRegions();

        forAll (remoteDonors, rdI)
        {
            if (remoteDonors[rdI] == index())
            {
                // Remote region lists current region as donor
                // Therefore, it is acceptor of this region
                accRegions[nAccRegions] = regionI;
                nAccRegions++;
            }
        }
    }

    // Reset size of acceptor regions
    accRegions.setSize(nAccRegions);
}


void Foam::oversetRegion::calcDonorAcceptorCells() const
{
    if (donorCellsPtr_ || acceptorCellsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcDonorAcceptorCells() const")
            << "Donor/acceptor cells already calculated"
            << abort(FatalError);
    }

    // If the pointers are not allocated, this means that the fringe has not
    // been calculated yet. Note that we need to perform fringe assembly for
    // each region in an iterative fashion because each layer of possibly new
    // acceptors will possibly invalidate eligible donors for this region.

    // Algorithm:
    //    ----- START ITERATIVE FRINGE ASSEMBLY PROCESS-----
    //    1) Loop through all regions and find donor/acceptor pairs for current
    //       set of acceptors as defined by oversetFringe,
    //    2) Check whether all regions have satisfied the Donor Suitability
    //       Criterion as defined by its oversetFringe
    //    3) If they have, stop the iterative process, otherwise repeat.
    //    4) Finalise overset assembly process by setting up donorCellsPtr_ and
    //       acceptorCellsPtr_ for all regions

    // Get all overset regions
    const PtrList<oversetRegion>& regions = oversetMesh_.regions();

    // Flag indicating that a suitable overlap has been found for all regions
    bool foundGlobalOverlap;

    do
    {
        // Set flag to true
        foundGlobalOverlap = true;

        // Loop through all regions
        forAll (regions, orI)
        {
            // Get current region
            const oversetRegion& curRegion = regions[orI];

            Info<< "Updating iteration for overset region: "
                << curRegion.name() << endl;

            // Update donor/acceptors for this region.
            // Note: updateDonorAcceptors() returns a bool indicating whether
            // a suitable overlap is found for this particular region
            const bool regionFoundSuitableOverlap =
                curRegion.updateDonorAcceptors();

            // Update global flag
            foundGlobalOverlap &= regionFoundSuitableOverlap;

            // If the overlap has not been found for this region, we need to
            // reset:
            //  - holeCells (depend on fringe holes)
            //  - eligibleDonors (depend on fringe holes and acceptors),
            //  - cellSearch (depends on eligible donors).
            if (!regionFoundSuitableOverlap)
            {
                deleteDemandDrivenData(curRegion.holeCellsPtr_);
                deleteDemandDrivenData(curRegion.eligibleDonorCellsPtr_);
                deleteDemandDrivenData(curRegion.cellSearchPtr_);
            }
        }
    } while (!foundGlobalOverlap);

    // Since a suitable overlap has been found for all regions, set-up
    // donor/acceptor fields for all regions
    forAll (regions, orI)
    {
        // Calling finaliseDonorAcceptor will take the latest set of suitable
        // donors/acceptors from its fringe algorithm and combine them into
        // donorCellsPtr_ and acceptorCellsPtr_
        regions[orI].finaliseDonorAcceptors();
    }
}


void Foam::oversetRegion::calcCutHoleCells() const
{
    if (cutHoleCellsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcCutHoleCells() const")
            << "Cut hole cells already calculated"
            << abort(FatalError);
    }

    // Algorithm
    // - go through all regions apart from the current region
    // - get access to hole surface search
    // - mark as hole all cells that fall outside of hole search
    // - combine all searches into a single inside-outside list

    // Get local cell indices
    const labelList& rc = regionCells();

    // Prepare local cell centres for inside-outside search on all regions
    vectorField localC(rc.size());

    const vectorField& c = mesh().cellCentres();

    forAll (localC, i)
    {
        localC[i] = c[rc[i]];
    }

    // Prepare hole mask
    boolList holeMask(mesh().nCells(), false);

    // Mark all hole cells using their hole boundary patch inside search

    // Get regions
    const PtrList<oversetRegion>& regions = oversetMesh_.regions();

    // Go through all regions apart from the current
    forAll (regions, regionI)
    {
        // Skip current region
        if (regionI == index())
        {
            continue;
        }

        const oversetRegion& otherRegion = regions[regionI];

        // If there are no hole patches on other region, skip it
        if (!otherRegion.holePatchesPresent())
        {
            continue;
        }

        // Get reference to hole search
        const triSurfaceSearch& holeSearch = otherRegion.holeSearch();

        boolList regionInside = holeSearch.calcInside(localC);

        // Note: hole mask has the size of all mesh cells and regionInside
        // only of the size of local region
        forAll (regionInside, i)
        {
            holeMask[rc[i]] |= regionInside[i];
        }
    }

    // Count hole cells
    label nHoleCells = 0;

    forAll (rc, i)
    {
        if (holeMask[rc[i]])
        {
            nHoleCells++;
        }
    }

    // Allocate hole cells storage
    cutHoleCellsPtr_ = new labelList(nHoleCells);
    labelList& ch = *cutHoleCellsPtr_;

    // Reset counter and collect hole cells
    nHoleCells = 0;

    forAll (rc, i)
    {
        if (holeMask[rc[i]])
        {
            ch[nHoleCells] = rc[i];
            nHoleCells++;
        }
    }

    if (oversetMesh::debug)
    {
        Pout<< "Region " << name()
            << " number of local holes = " << cutHoleCellsPtr_->size()
            << endl;
    }
}


void Foam::oversetRegion::calcHoleCells() const
{
    if (holeCellsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcHoleCells() const")
            << "Hole cells already calculated"
            << abort(FatalError);
    }

    // Combine cut hole cells and fringe cells into a single list

    // Prepare hole mask
    boolList holeMask(mesh().nCells(), false);

    // Mask all cut hole cells
    const labelList& cutHoleCells = cutHoles();

    forAll (cutHoleCells, i)
    {
        holeMask[cutHoleCells[i]] = true;
    }

    // Mask fringe hole cells
    const labelList& fringeHoleCells = fringePtr_->fringeHoles();

    forAll (fringeHoleCells, i)
    {
        holeMask[fringeHoleCells[i]] = true;
    }

    // Count hole cells in the region
    const labelList& rc = regionCells();

    label nHoleCells = 0;

    forAll (rc, i)
    {
        if (holeMask[rc[i]])
        {
            ++nHoleCells;
        }
    }

    // Allocate hole cells storage
    holeCellsPtr_ = new labelList(nHoleCells);
    labelList& h = *holeCellsPtr_;

    // Reset counter and collect hole cells
    nHoleCells = 0;

    forAll (rc, i)
    {
        if (holeMask[rc[i]])
        {
            h[nHoleCells] = rc[i];
            ++nHoleCells;
        }
    }
}


void Foam::oversetRegion::calcEligibleDonorCells() const
{
    if (eligibleDonorCellsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcEligibleDonorCells() const")
            << "Eligible donor cells already calculated"
            << abort(FatalError);
    }

    const labelList& rc = regionCells();

    // Prepare donor mask for region cells
    boolList donorMask(mesh().nCells(), false);

    // Mark region cells as eligible
    forAll (rc, rcI)
    {
        donorMask[rc[rcI]] = true;
    }

    // Remove all hole cells from region
    const labelList& hc = holes();

    forAll (hc, hcI)
    {
        donorMask[hc[hcI]] = false;
    }

    // Remove all acceptor cells from fringe
    const labelList& ac = fringePtr_->candidateAcceptors();

    forAll (ac, acI)
    {
        donorMask[ac[acI]] = false;
    }

    // Count and collect eligible donors
    label nEligibleDonors = 0;

    forAll (donorMask, cellI)
    {
        if (donorMask[cellI])
        {
            nEligibleDonors++;
        }
    }

    eligibleDonorCellsPtr_ = new labelList(nEligibleDonors);
    labelList& ed = *eligibleDonorCellsPtr_;

    // Reset counter and collect cells
    nEligibleDonors = 0;

    forAll (donorMask, cellI)
    {
        if (donorMask[cellI])
        {
            ed[nEligibleDonors] = cellI;
            nEligibleDonors++;
        }
    }
}


void Foam::oversetRegion::calcHoleTriMesh() const
{
    if (holeTriMeshPtr_)
    {
        FatalErrorIn("void oversetRegion::calcHoleTriMesh() const")
            << "Hole tri mesh already calculated"
            << abort(FatalError);
    }

    // Create region mask to check if patch touches region
    boolList regionMask(mesh().nCells(), false);

    const labelList& rc = regionCells();

    forAll (rc, rcI)
    {
        regionMask[rc[rcI]] = true;
    }

    // Get hole patch names
    const wordList& holePatchNames = overset().holePatchNames();

    // Collect local hole faces
    labelHashSet holePatches;

    forAll (holePatchNames, nameI)
    {
        polyPatchID curHolePatch
        (
            holePatchNames[nameI],
            mesh().boundaryMesh()
        );

        if (curHolePatch.active())
        {
            // If the patch has zero size, do not insert it
            // Parallel cutting bug.  HJ, 17/Apr/2014
            if (!mesh().boundaryMesh()[curHolePatch.index()].empty())
            {
                // Check if the patch is touching the current region
                const labelList& faceCells =
                    mesh().boundary()[curHolePatch.index()].faceCells();

                label nFound = 0;

                forAll (faceCells, fcI)
                {
                    if (regionMask[faceCells[fcI]])
                    {
                        nFound++;
                    }
                }

                // Check if the complete patch belongs to current region
                if (nFound == faceCells.size())
                {
                    holePatches.insert(curHolePatch.index());
                }
                else if (nFound > 0)
                {
                    WarningIn("void oversetRegion::calcHoleTriMesh() const")
                        << "Patch " << holePatchNames[nameI]
                        << " seems to be split between multiple regions.  "
                        << "Please check overset region structure.  "
                        << "nFound: " << nFound
                        << " faceCells: " << faceCells.size()
                        << endl;
                }
            }
        }
        else
        {
            FatalErrorIn
            (
                "const triSurfaceMesh& oversetRegion::holeTriMesh() const"
            )   << "Patch "  << holePatchNames[nameI]
                << " cannot be found.  Available patch names: "
                << mesh().boundaryMesh().names()
                << abort(FatalError);
        }
    }

    // Make and invert local triSurface
    triFaceList triFaces;
    pointField triPoints;

    // Memory management
    {
        triSurface ts = triSurfaceTools::triangulate
        (
            mesh().boundaryMesh(),
            holePatches
        );

        // Clean mutiple points and zero-sized triangles
        ts.cleanup(false);

        triFaces.setSize(ts.size());
        triPoints = ts.points();

        forAll (ts, tsI)
        {
            triFaces[tsI] = ts[tsI].reverseFace();
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
                    renumberFace[rfI] = procRenumberPoints[procFace[rfI]];
                }

                nTris++;
            }
        }
    }

    // Make a complete triSurface from local data
    holeTriMeshPtr_ = new triSurface
    (
        triFaces,
        triPoints
    );

    // Clean up duplicate points and zero sized triangles
    holeTriMeshPtr_->cleanup(false);

    Info<< "Region " << name() << ": "
        << holeTriMeshPtr_->size() << " triangles in hole cutting"
        << endl;

    // Debug: write holeTriMesh
    if (Pstream::master())
    {
        if (!holeTriMeshPtr_->empty())
        {
            holeTriMeshPtr_->write(word("holeTriSurface_") + name() + ".vtk");
        }
    }
}


void Foam::oversetRegion::calcBounds() const
{
    if (localBoundsPtr_ || globalBoundsPtr_)
    {
        FatalErrorIn("void oversetRegion::calcBounds() const")
            << "Bounds already calculated"
            << abort(FatalError);
    }

    // Make a global bounding box for this region
    boolList usedPoints(mesh_.nPoints(), false);

    // Get cells-points from mesh
    const labelListList& pc = mesh_.cellPoints();

    // Get region cell indices
    const labelList& rc = zone();

    forAll (rc, rcI)
    {
        // Get points of region cells
        const labelList& curPc = pc[rc[rcI]];

        forAll (curPc, i)
        {
            usedPoints[curPc[i]] = true;
        }
    }

    // Count used points
    label nUsedPoints = 0;

    forAll (usedPoints, pointI)
    {
        if (usedPoints[pointI])
        {
            nUsedPoints++;
        }
    }

    // Make a list of used points
    const pointField& points = mesh_.points();

    pointField regionPoints(nUsedPoints);

    // Reset point counter to zero
    nUsedPoints = 0;

    forAll (usedPoints, pointI)
    {
        if (usedPoints[pointI])
        {
            regionPoints[nUsedPoints] = points[pointI];
            nUsedPoints++;
        }
    }

    // Local (processor) bounding box is calculated without a reduce
    localBoundsPtr_ = new boundBox(regionPoints, false);

    // Global bounding box is calculated with a reduce
    globalBoundsPtr_ = new boundBox(regionPoints, true);
}


void Foam::oversetRegion::calcCellSearch() const
{
    if (cellSearchPtr_)
    {
        FatalErrorIn("void oversetRegion::calcCellSearch() const")
            << "Cell tree already calculated"
            << abort(FatalError);
    }

    // Create the octree search for this region.  It will be used by other
    // regions when searching for donor cells

    // Bounding box containing only local region cells
    treeBoundBox overallBb(localBounds());
    Random rndGen(123456);
    overallBb = overallBb.extend(rndGen, 1E-4);
    overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    // Search
    cellSearchPtr_ = new indexedOctree<treeDataCell>
    (
        treeDataCell
        (
            false,  //  Cache bb.  Reconsider for moving mesh cases
            mesh_,
            eligibleDonors()
        ),
        overallBb,  // overall search domain
        8,          // maxLevel
        10,         // leafsize
        3.0         // duplicity
    );
}


void Foam::oversetRegion::calcProcBoundBoxes() const
{
    if (procBoundBoxesPtr_)
    {
        FatalErrorIn("void oversetRegion::calcProcBoundBoxes() const")
            << "Processor bounding boxes already calculated"
            << abort(FatalError);
    }

    // Create the list
    procBoundBoxesPtr_ = new List<List<boundBox> >(Pstream::nProcs());
    List<List<boundBox> >& procBoundBoxes = *procBoundBoxesPtr_;

    // Get pointer list of bounding boxes for this processor
    List<boundBox>& localBoundBoxes = procBoundBoxes[Pstream::myProcNo()];

    // Get all regions that are present on this processor
    const PtrList<oversetRegion>& regions = oversetMesh_.regions();

    // Set the size for this processor
    localBoundBoxes.setSize(regions.size());

    // Loop through overset regions and populate the list
    forAll (regions, orI)
    {
        if (useLocalBoundBoxes_)
        {
            // Use local bounds to optimise sending of acceptors
            localBoundBoxes[orI] = regions[orI].localBounds();
        }
        else
        {
            // Use global bounds to make sure we find donors to all acceptors,
            // even if the bounding boxes of certain regions do not overlap
            localBoundBoxes[orI] = this->globalBounds();
        }
    }

    // Now that each processor has filled in its own part, combine the data
    Pstream::gatherList(procBoundBoxes);
    Pstream::scatterList(procBoundBoxes);
}


void Foam::oversetRegion::clearOut() const
{
    deleteDemandDrivenData(donorRegionsPtr_);
    deleteDemandDrivenData(acceptorRegionsPtr_);

    deleteDemandDrivenData(acceptorCellsPtr_);
    deleteDemandDrivenData(donorCellsPtr_);
    deleteDemandDrivenData(cutHoleCellsPtr_);
    deleteDemandDrivenData(holeCellsPtr_);
    deleteDemandDrivenData(eligibleDonorCellsPtr_);

    deleteDemandDrivenData(holeTriMeshPtr_);
    deleteDemandDrivenData(holeSearchPtr_);

    deleteDemandDrivenData(localBoundsPtr_);
    deleteDemandDrivenData(globalBoundsPtr_);
    deleteDemandDrivenData(cellSearchPtr_);

    deleteDemandDrivenData(procBoundBoxesPtr_);
}


bool Foam::oversetRegion::updateDonorAcceptors() const
{
    // If a suitable fringe on this region has been found, simply return true
    if (fringePtr_->foundSuitableOverlap())
    {
        return true;
    }

    // Algorithm (operates on acceptor cells):
    // 1) Start iteration for this region
    //
    //    -----SENDING/RECEIVING ACCEPTORS PART-----
    //
    // 2) We need to calculate the sending map for acceptors, which tells me
    //    which acceptors I need to send to which processor:
    //    - Loop through acceptors for this region and through all donor regions
    //    - Using the processor bounding boxes, figure out where to send this
    //      acceptor (there is no need to send the acceptor if the acceptor
    //      point falls outside the donor region of a certain processor). For
    //      each processor, simply insert the acceptor index that needs sending
    //      into the DynamicList
    //    - While looping through acceptors, count how many acceptors I'm
    //      sending to each processor
    // 3) Count how many acceptors my (local) processor is going to be
    //    receiving from all other processors
    // 4) Create a constructing map for acceptors, organized as follows:
    //    - If processor N sends me M acceptor points, these points will be
    //      stored in the receiving list from indices K to K + M - 1, where K
    //      is the total number of acceptor points received from previous
    //      processors (processors I, where I < N)
    // 5) Create mapDistribute object and distribute acceptor data
    //
    //    -----DONOR SEARCH PART-----
    //
    // 6) Find donors for received acceptors (using octree search based on
    //    eligible donors only) and collect them in donorAcceptor list
    //    necessary for addressing after communicating donor data
    //
    //    -----SENDING/RECEIVING DONORS PART-----
    //
    // 7) Create sending map for donors - this is simply the constructing map
    //    for acceptors
    // 8) Count how many donor/acceptor pairs my (local) processor is going to
    //    be receiving from all other processors. Note that for a given
    //    acceptor we can actually have multiple donors coming from different
    //    processors
    // 9) Create constructing map for donors, organized as follows (the same way
    //    as the constructing map for acceptors):
    //    - If processor N sends me M donors, these donors will be stored in the
    //      receiving list from indices K to K + M - 1, where K is the total
    //      number of donor points received from previous processors (processors
    //      I, where I < N)
    // 10) Create mapDistribute object and distribute donor data
    // 11) Filter possibly multiple remote donors (coming from different
    //     processors) for all acceptors
    //
    //    -----FRINGE HANDLING-----
    //
    // 12) Pass new batch of donor/acceptor pairs to fringe algorithm,
    //     which does its own magic (donor suitability, iteration control,
    //     update of eligible donors...)

    // Get necessary mesh data
    const vectorField& cc = mesh_.cellCentres();
    const labelListList& cCells = mesh_.cellCells(); // For extended donors

    // Get donor regions
    const labelList& dr = donorRegions();

    // Get all overset regions
    const PtrList<oversetRegion>& regions = oversetMesh_.regions();

    // Get region bounding boxes for all processors and regions: used to decide
    // where we need to send the data (i.e. it does not make sense to send
    // acceptor data to a processor if the acceptor point does not lie within
    // the bounding box of a given processor donor region)
    const List<List<boundBox> >& procRegionBB = procBoundBoxes();

    // STAGE 1: Start iteration for this region

    // Get current and local acceptor cells from oversetFringe
    const labelList& a = fringePtr_->candidateAcceptors();

    // Create a list of local acceptors (holding acceptors on this processor and
    // donors on possibly remote processor)
    donorAcceptorList localAcceptorDonorList(a.size());

    // Insert local acceptor into the list
    forAll (localAcceptorDonorList, aI)
    {
        localAcceptorDonorList[aI] = donorAcceptor
        (
            aI, // Note: storing index in the local acceptor list (not the
                // cell index!). Done this way for easier filtering of
                // multiple donors in STAGE 11
            Pstream::myProcNo(),
            cc[a[aI]]
        );
    }

    // Create list containing number of acceptors my processor is sending to
    // other processors.
    // Example: nAcceptorsToProcessorMap[2][5] = 123 will (after collection
    // of data) mean that processor 2 sends 123 acceptor points to processor
    // 5. Note that for a serial run, this list is completely unnecessary,
    // but I prefer writing this in a general way, where I don't care about
    // minor loss of efficiency for serial runs. VV, 30/Jan/2016.
    labelListList nAcceptorsToProcessorMap(Pstream::nProcs());

    forAll (nAcceptorsToProcessorMap, procI)
    {
        nAcceptorsToProcessorMap[procI].setSize(Pstream::nProcs(), 0);
    }

    // Get number of acceptors I'm sending to other processors
    labelList& numberOfLocalAcceptorsToProcs =
        nAcceptorsToProcessorMap[Pstream::myProcNo()];

    // STAGE 2: Calculate the sending map (for sending acceptor points)

    // Example: sendAcceptorMap[procI] = (0 1 5 7 89 ...) tells us that I
    // should send field values (acceptor points in this case) indexed by:
    // 0, 1, 5, 7, 89... to processor procI

    // Initialize sending map: for each processor, create a DynamicList of
    // acceptors that need to be sent to that processor
    List<dynamicLabelList> sendAcceptorMap(Pstream::nProcs());

    // Allocate enough storage as if we are sending all acceptors to all
    // processors (trading off memory for performance)
    forAll (sendAcceptorMap, procI)
    {
        sendAcceptorMap[procI].setCapacity(a.size());
    }

    // Loop through all processors
    forAll (sendAcceptorMap, procI)
    {
        // Get bounding boxes on this processor
        const List<boundBox>& curProcBoundBoxes =
            procRegionBB[procI];

        // Get current processor send map
        dynamicLabelList& curSendMap = sendAcceptorMap[procI];

        // Loop through all donor regions
        forAll (dr, drI)
        {
            // Get region index of this donor region
            const label& curDonorRegion = dr[drI];

            // Loop through all local acceptors
            forAll (a, aI)
            {
                // Check whether the acceptor is within the bounding box of
                // this donor region on this particular processor.
                if
                (
                    curProcBoundBoxes[curDonorRegion].containsInside
                    (
                        localAcceptorDonorList[aI].acceptorPoint()
                    )
                )
                {
                    // Acceptor may find donor on this processor, append it
                    curSendMap.append(aI);

                    // Increment the number of acceptors I'm sending to this
                    // processor
                    ++numberOfLocalAcceptorsToProcs[procI];
                }
            } // End for all local acceptors
        } // End for all donor regions
    } // End for all processors

    // STAGE 3: Count number of points I'm receiving from all other
    // processors

    // Gather/scatter number of acceptor points going to each processor from
    // each processor so that all processors have all necessary information
    // when creating the map distribute tool for distributing acceptor
    // points
    Pstream::gatherList(nAcceptorsToProcessorMap);
    Pstream::scatterList(nAcceptorsToProcessorMap);

    // Count how many acceptors I'm going to receive from others
    label nAcceptorReceives = 0;
    forAll (nAcceptorsToProcessorMap, procI)
    {
        nAcceptorReceives +=
            nAcceptorsToProcessorMap[procI][Pstream::myProcNo()];
    }

    // STAGE 4: Calculation of construct map for acceptors

    // The construct map is simply an index offseted by the number of values
    // received by previous processors.
    // Example:
    /*
        Procs sending to me | Number of items being sent
        ------------------------------------------------
               P0           |             1
               P1           |             7
               P5           |             2
               .            |             .
               .            |             .
               .            |             .

        Received data has the following form:
        (
            a_0, (one value from proc 0)
            a_1, a_2, a_3, a_4, a_5, a_6, a_7 (seven values from proc 1)
            a_8, a_9 (two values from proc 5)
            ...
            ...
            ...
        )
    */

    // Create construct map
    labelListList constructAcceptorMap(Pstream::nProcs());

    // Counter for offset
    label procOffset = 0;

    forAll (constructAcceptorMap, procI)
    {
        // Get receiving size from this processor
        const label nReceivesFromCurProc =
            nAcceptorsToProcessorMap[procI][Pstream::myProcNo()];

        // Get current construct map
        labelList& curConstructMap = constructAcceptorMap[procI];

        // Set the size corresponding to number of received acceptor points
        curConstructMap.setSize(nReceivesFromCurProc);

        // Set mapping as a simple offset
        forAll (curConstructMap, receivedItemI)
        {
            curConstructMap[receivedItemI] = receivedItemI + procOffset;
        }

        // Increment the processor offset by the size received from this
        // processor
        procOffset += nReceivesFromCurProc;
    }

    // STAGE 5: Distribute acceptor points

    // Need to create a labelListList from List<dynamicLabelList> for sending
    // map.
    labelListList sendAcceptorFixedMap(Pstream::nProcs());
    forAll (sendAcceptorFixedMap, procI)
    {
        // Transfer the content of dynamic list into this processor list,
        // sendAcceptorMap is invalid from now on
        sendAcceptorFixedMap[procI].transfer(sendAcceptorMap[procI]);
    }

    // Create mapDistribute object for distributing acceptor points. Note:
    // reusing maps, meaning that arguments are invalid from now onward.
    mapDistribute acceptorDistribution
    (
        nAcceptorReceives,
        sendAcceptorFixedMap,
        constructAcceptorMap,
        true // reuse maps
    );

    // Distribute acceptor data. Note: now localAcceptorDonorList holds
    // acceptor data received from other processor for which we need to find
    // eligible donors
    acceptorDistribution.distribute(localAcceptorDonorList);

    // Use an alias (reference) from now on for clarity
    donorAcceptorList& receivedAcceptorDonorList = localAcceptorDonorList;

    // STAGE 6: Find donors for received acceptors

    // For each donor region, create a list of donor/acceptor pairs
    // data.
    // Note 1: the addressing of donor/acceptors will be the same as
    // received acceptors in order to easily send the data back in case of
    // multiple donor regions.
    // Note 2: We will prefer donors that are closer to the acceptor and
    // send back only those.

    // Loop through donor regions
    forAll (dr, drI)
    {
        // Get current donor region
        const oversetRegion& curDonorRegion = regions[dr[drI]];

        // Get current eligible donors for this region. Note: these are
        // updated after a particular overset region finishes an iteration
        // of fringe assembly
        const labelList& curDonors = curDonorRegion.eligibleDonors();

        // Mask eligible donors for extended neighbourhood search
        boolList eligibleDonorMask(mesh_.nCells(), false);
        forAll (curDonors, i)
        {
            eligibleDonorMask[curDonors[i]] = true;
        }

        // Get donor region tree. Note: octree also depends on eligible
        // donors, updated at the end of iteration
        const indexedOctree<treeDataCell>& tree =
            curDonorRegion.cellSearch();

        // It is possible that an octree is empty (if there are no eligible
        // donor cells on this processor), do not search
        if (tree.nodes().empty())
        {
            continue;
        }

        const scalar span = tree.bb().mag();

        // Loop through received acceptor data
        forAll (receivedAcceptorDonorList, accI)
        {
            // Get current donor/acceptor pair
            donorAcceptor& daPair = receivedAcceptorDonorList[accI];

            // Get acceptor cell centre
            const point& curP = daPair.acceptorPoint();

            // Find nearest donor cell with octree. Note: octree only
            // contains eligible cells.  HJ, 10/Jan/2015.
            const pointIndexHit pih = tree.findNearest(curP, span);

            if (pih.hit())
            {
                // Found a hit, check whether this donor is set or not. If
                // it is not set, set it no questions asked; if it is set,
                // check whether this is a better candidate by either
                // looking whether acceptor point is within donor cell or
                // taking a closer hit

                // Get index obtained by octree
                const label donorCandidateIndex = pih.index();

                if
                (
                   !daPair.donorFound()
                 || mesh_.pointInCellBB
                    (
                        curP,
                        curDonors[donorCandidateIndex]
                    )
                 || (
                        mag(cc[curDonors[donorCandidateIndex]] - curP)
                      < mag(daPair.donorPoint() - curP)
                    )
                )
                {
                    // Set donor
                    daPair.setDonor
                    (
                        curDonors[donorCandidateIndex],
                        Pstream::myProcNo(),
                        cc[curDonors[donorCandidateIndex]]
                    );

                    // Set extended donors
                    daPair.setExtendedDonors
                    (
                        eligibleDonorMask,
                        cCells,
                        cc
                    );
                }

                // Note: consider removing pointInCellBB since it now has
                // precedence over distance criterion.  VV, 31/Jan/2017.
            }
            else if (oversetMesh::debug && !daPair.donorFound())
            {
                // This donor is not valid and I did not find a hit in
                // octree, issue a warning
                WarningIn
                (
                    "void oversetRegion::updateDonorAcceptors() const"
                )   << "Could not find a hit for acceptor,"
                    << "donor may remain invalid."
                    << endl;
            }
        } // End for all acceptor cell centres
    } // End for all donor regions

    // STAGE 7: Create sending map for donors

    // Note: sending map for donors is basically the constructing map for
    // acceptors since we have used the same addressing for donor search.
    // Create a copy from map distribute object used to communicate acceptor
    // data
    labelListList sendDonorMap = acceptorDistribution.constructMap();

    // STAGE 8: Count how many donors I'm going to receive

    // Note: reuse nAcceptorsToProcessorMap[i][j], which tells me how many
    // acceptors processor i is sending to processor j. The number of donors
    // received by processor j is the same as the number of acceptors sent
    // to processor j.
    label nDonorReceives = 0;
    forAll(nAcceptorsToProcessorMap, procI)
    {
        nDonorReceives +=
            nAcceptorsToProcessorMap[Pstream::myProcNo()][procI];
    }

    // STAGE 9: Calculation of construct map for donors

    // Create construct map
    labelListList constructDonorMap(Pstream::nProcs());

    // Reset existing processor offset counter
    procOffset = 0;

    forAll (constructDonorMap, procI)
    {
        // Get receiving size from this processor. Reusing
        // nAcceptorsToProcessorMap as in STAGE 8
        const label nReceivesFromCurProc =
            nAcceptorsToProcessorMap[Pstream::myProcNo()][procI];

        // Get current construct map
        labelList& curConstructMap = constructDonorMap[procI];

        // Set the size corresponding to number of received donor/acceptor
        // pairs
        curConstructMap.setSize(nReceivesFromCurProc);

        // Set mapping as a simple offset
        forAll (curConstructMap, receivedItemI)
        {
            curConstructMap[receivedItemI] = receivedItemI + procOffset;
        }

        // Increment the processor offset by the size received from this
        // processor
        procOffset += nReceivesFromCurProc;
    }

    // STAGE 10: Distribute donor data

    // Create mapDistribute object for distributing donor data. Note:
    // reusing maps, meaning that arguments are invalid from now onward.
    mapDistribute donorDistribution
    (
        nDonorReceives,
        sendDonorMap,
        constructDonorMap,
        true // reuse maps
    );

    // Distribute donor/acceptor pairs, for initial N acceptors I have sent
    // across certain processors, I will receive M donor/acceptor pairs,
    // where M >= N.
    donorDistribution.distribute(receivedAcceptorDonorList);

    // Use an alias (reference) from now on for clarity
    donorAcceptorList& completeDonorAcceptorList = receivedAcceptorDonorList;

    // Before filtering, check whether all received donors are actually for
    // acceptors on this processor. If not, something went terribly wrong. Used
    // for testing/debugging parallel comms
    if (oversetMesh::debug)
    {
        forAll (completeDonorAcceptorList, daI)
        {
            if
            (
                completeDonorAcceptorList[daI].acceptorProcNo()
             != Pstream::myProcNo()
            )
            {
                FatalErrorIn("void oversetRegion::updateDonorAcceptors() const")
                    << "Received donor/acceptor pair where acceptor belongs to "
                    << "a different processor. " << nl
                    << "My processor number: " << Pstream::myProcNo()
                    << "Acceptor processor number: "
                    << completeDonorAcceptorList[daI].acceptorProcNo()
                    << abort(FatalError);
            }
        }
    }

    // STAGE 11: Filter possibly multiple remote donors

    // Create a masking field indicating that a certain acceptor has been
    // visited
    boolList isVisited(a.size(), false);

    // Create a combined donor acceptor list, only containing best donors for
    // current acceptors.
    donorAcceptorList combinedDonorAcceptorList(a.size());

    // Loop through donor/acceptor list
    forAll (completeDonorAcceptorList, daI)
    {
        // Get current donor/acceptor pair
        const donorAcceptor& curDA = completeDonorAcceptorList[daI];

        // Get current acceptor index (not the cell index, but the index into
        // the current acceptor list). See STAGE 1
        const label& aI = curDA.acceptorCell();

        // Get the flag indicating whether the acceptor has been visited or not
        bool& acceptorVisited = isVisited[aI];

        // Get current combined donor/acceptor pair
        donorAcceptor& curDACombined = combinedDonorAcceptorList[aI];

        if (!acceptorVisited)
        {
            // This acceptor has not been previously visited, set it in the
            // combined list
            curDACombined = curDA;

            // Set the correct cell index in the combined list
            curDACombined.acceptorCell() = a[aI];

            // Mark as visited
            acceptorVisited = true;
        }
        else
        {
            // This acceptor has been previously visited, meaning we have to
            // make a choice whether to update it or not. At this point, the
            // choice will be based on least distance from acceptor cell centre
            // to donor cell centre. Run-time selectable Donor Suitability
            // Functions will be applied in oversetFringe
            if (curDA.distance() < curDACombined.distance())
            {
                // This is a better candidate for the same acceptor, set donor
                // accordingly
                curDACombined.setDonor
                (
                    curDA.donorCell(),
                    curDA.donorProcNo(),
                    curDA.donorPoint()
                );
            }
        }
    }

    // Check whether all acceptors have been visited. Used for testing/debugging
    // parallel comms
    if (oversetMesh::debug)
    {
        bool allVisited = true;

        forAll (isVisited, aI)
        {
            allVisited &= isVisited[aI];
        }

        if (!allVisited)
        {
            if (!useLocalBoundBoxes_)
            {
                FatalErrorIn("void oversetRegion::updateDonorAcceptors() const")
                    << "Did not visit all acceptors when recombining data..."
                    << nl
                    << "... and we did not use local processor bounding boxes."
                    << nl
                    << "Something went wrong."
                    << abort(FatalError);
            }
            else
            {
                FatalErrorIn("void oversetRegion::updateDonorAcceptors() const")
                    << "Did not visit all acceptors when recombining data..."
                    << nl
                    << "Try switching off useLocalBoundingBoxes for all regions"
                    << nl
                    << "(this optimisation is switched on by default)."
                    << abort(FatalError);
            }
        }
    }

    // STAGE 12: Finish the iteration by updating the fringe, which will
    // actually hold final and some intermediate steps for donor/acceptor
    // assembly process

    // Note, may invalidate the argument list depending on the fringe
    // algorithm that is used
    bool suitableOverlapFound =
        fringePtr_->updateIteration(combinedDonorAcceptorList);

    return suitableOverlapFound;
}


void Foam::oversetRegion::finaliseDonorAcceptors() const
{
    // Check the state of pointers (if someone calls this function from some
    // member function other than calcDonorAcceptorCells)
    if (donorCellsPtr_ || acceptorCellsPtr_)
    {
        FatalErrorIn("void oversetRegion::finaliseDonorAcceptors() const")
            << "Donor/acceptor cells already calculated. Make sure you have "
            << "called this function from calcAcceptorsDonors only."
            << abort(FatalError);
    }

    // Need to fetch final donor/acceptor pair from fringe and recombine donor
    // cells for parallel run (donorCellsPtr_ will hold a list containing cells
    // where the donor is local and acceptor is possibly remote)

    // Algorithm:
    // 1) Acceptor cells are simply donor/acceptor pairs from final fringe
    //    iteration. Get them.
    // 2) Calculate the sending map, telling me how many donor/acceptor pairs I
    //    need to send to a certain processor. The algorithm is essentially the
    //    same as STAGE 2 in oversetRegion::updateDonorAcceptors()
    // 3) Count how many donor/acceptor pairs my (local) processor is going to
    //    be receiving from all other processors
    // 4) Create a constructing map for donor/acceptor pairs, organized as in
    //    STAGE 4 in oversetRegion::updateDonorAcceptors()
    // 5) Create mapDistribute object and distribute local acceptor data, thus
    //    creating local donor data

    // STAGE 1: Get acceptor cells

    // Reuse the list from fringe handler (thus invalidating it)
    acceptorCellsPtr_ = new donorAcceptorList
    (
        fringePtr_->finalDonorAcceptors(),
        true // reuse
    );
    const donorAcceptorList& acceptorCells = *acceptorCellsPtr_;

    // Note: the only way it is possible that the donor has not been found is if
    // the original acceptor cell is not within bounding boxes of all processor
    // containing donor regions. This means that the mesh set-up is considered
    // invalid, although this will not fail on serial runs because we do not
    // require that the acceptor cell is within the donor cell. VV, 6/Feb/2017.
    forAll (acceptorCells, accI)
    {
        if (!acceptorCells[accI].donorFound())
        {
            FatalErrorIn("void oversetRegion::finaliseDonorAcceptors() const")
                << "Did not find a donor for acceptor at: "
                << acceptorCells[accI].acceptorPoint()
                << nl
                << "This means that all donor regions do not contain acceptor "
                << "point, implying invalid overset mesh."
                << nl
                << "Please check your overset mesh structure."
                << abort(FatalError);
        }
    }

    // STAGE 2: Calculate the sending map

    // Create list containing number of donor/acceptor pairs my processor is
    // sending to other processors
    labelListList nPairsToProcessorMap(Pstream::nProcs());

    forAll (nPairsToProcessorMap, procI)
    {
        nPairsToProcessorMap[procI].setSize(Pstream::nProcs(), 0);
    }

    // Get local list (number of acceptors I'm sending to other processors)
    labelList& numberOfLocalAcceptorPairsToProcs =
        nPairsToProcessorMap[Pstream::myProcNo()];

    // Initialize sending map: for each processor, create a DynamicList of
    // donor/acceptor pairs that need to be sent to that processor
    List<dynamicLabelList> sendMap(Pstream::nProcs());

    // Allocate enough storage as if we are sending all pairs to all
    // processors (trading off memory for performance)
    forAll (sendMap, procI)
    {
        sendMap[procI].setCapacity(acceptorCells.size());
    }

    // Loop through all donor/acceptor pairs where acceptor is on the local
    // processor
    forAll (acceptorCells, aI)
    {
        // Get donor processor index I need to send data to
        const label& donorProcID = acceptorCells[aI].donorProcNo();

        // Append index to the sending map for donor processor
        sendMap[donorProcID].append(aI);

        // Increment the number of acceptors I'm sending to this processor
        ++numberOfLocalAcceptorPairsToProcs[donorProcID];
    }

    // STAGE 3: Count number of donor/acceptor pairs I'm receiving from all
    // other processors

    // Gather/scatter in order to have complete data: how many donor/acceptor
    // pairs are actually sent from each processor to all other processors
    Pstream::gatherList(nPairsToProcessorMap);
    Pstream::scatterList(nPairsToProcessorMap);

    // Count how many donor/acceptor pairs I'm going to receive from others
    label nReceives = 0;
    forAll (nPairsToProcessorMap, procI)
    {
        nReceives += nPairsToProcessorMap[procI][Pstream::myProcNo()];
    }

    // STAGE 4: Calculation of construct map: index offseted by the number of
    // values received by previous processors. See STAGE 4 in
    // oversetRegion::updateDonorsAcceptors() for details

    // Create construct map
    labelListList constructMap(Pstream::nProcs());

    // Counter for offset
    label procOffset = 0;

    forAll (constructMap, procI)
    {
        // Get receiving size from this processor
        const label nReceivesFromCurProc =
            nPairsToProcessorMap[procI][Pstream::myProcNo()];

        // Get current construct map
        labelList& curConstructMap = constructMap[procI];

        // Set the size corresponding to number of received pairs
        curConstructMap.setSize(nReceivesFromCurProc);

        // Set mapping as a simple offset
        forAll (curConstructMap, receivedItemI)
        {
            curConstructMap[receivedItemI] = receivedItemI + procOffset;
        }

        // Increment the processor offset by the size received from this
        // processor
        procOffset += nReceivesFromCurProc;
    }

    // STAGE 5: Distribute donor/acceptor pairs

    // Need to create a labelListList from List<dynamicLabelList> for sending
    // map.
    labelListList sendFixedMap(Pstream::nProcs());
    forAll (sendFixedMap, procI)
    {
        // Transfer the content of dynamic list into this processor list,
        // sendMap is invalid from now on
        sendFixedMap[procI].transfer(sendMap[procI]);
    }

    // Create mapDistribute object. Note: reusing maps, meaning that arguments
    // are invalid from now onward.
    mapDistribute localAcceptorToLocalDonor
    (
        nReceives,
        sendFixedMap,
        constructMap,
        true // reuse maps
    );

    // Copy local acceptor data to initialize local donor data
    donorCellsPtr_ = new donorAcceptorList(acceptorCells);

    // Distribute local acceptor data, thus creating local donor data after
    // distribution
    localAcceptorToLocalDonor.distribute(*donorCellsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetRegion::oversetRegion
(
    const word& name,
    const label index,
    const fvMesh& mesh,
    const oversetMesh& oversetMesh,
    const dictionary& dict
)
:
    name_(name),
    index_(index),
    mesh_(mesh),
    oversetMesh_(oversetMesh),
    zoneIndex_(mesh_.cellZones().findZoneID(name_)),
    donorRegionNames_(dict.lookup("donorRegions")),
    fringePtr_(),
    donorRegionsPtr_(NULL),
    acceptorRegionsPtr_(NULL),

    acceptorCellsPtr_(NULL),
    donorCellsPtr_(NULL),
    cutHoleCellsPtr_(NULL),
    holeCellsPtr_(NULL),
    eligibleDonorCellsPtr_(NULL),

    holeTriMeshPtr_(NULL),
    holeSearchPtr_(NULL),

    localBoundsPtr_(NULL),
    globalBoundsPtr_(NULL),
    cellSearchPtr_(NULL),
    procBoundBoxesPtr_(NULL),
    useLocalBoundBoxes_
    (
        dict.lookupOrDefault<Switch>
        (
            "useLocalBoundBoxes",
            false
        )
    )
{
    // Check zone index
    if (zoneIndex_ < 0)
    {
        FatalErrorIn
        (
            "oversetRegion::oversetRegion\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const oversetMesh& oversetMesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Cannot find cell zone for region " << name << nl
            << "Available cell zones: " << mesh_.cellZones().names()
            << abort(FatalError);
    }

    fringePtr_ = oversetFringe::New
    (
        mesh,
        *this,
        dict.subDict("fringe")
    );

    calcBounds();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::oversetRegion::~oversetRegion()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::labelList& Foam::oversetRegion::donorRegions() const
{
    if (!donorRegionsPtr_)
    {
        calcDonorRegions();
    }

    return *donorRegionsPtr_;
}


const Foam::labelList& Foam::oversetRegion::acceptorRegions() const
{
    if (!acceptorRegionsPtr_)
    {
        calcAcceptorRegions();
    }

    return *acceptorRegionsPtr_;
}


const Foam::donorAcceptorList& Foam::oversetRegion::acceptors() const
{
    if (!acceptorCellsPtr_)
    {
        calcDonorAcceptorCells();
    }

    return *acceptorCellsPtr_;
}


const Foam::donorAcceptorList& Foam::oversetRegion::donors() const
{
    if (!donorCellsPtr_)
    {
        calcDonorAcceptorCells();
    }

    return *donorCellsPtr_;
}


const Foam::labelList& Foam::oversetRegion::cutHoles() const
{
    if (!cutHoleCellsPtr_)
    {
        calcCutHoleCells();
    }

    return *cutHoleCellsPtr_;
}


const Foam::labelList& Foam::oversetRegion::holes() const
{
    if (!holeCellsPtr_)
    {
        calcHoleCells();
    }

    return *holeCellsPtr_;
}


const Foam::labelList& Foam::oversetRegion::eligibleDonors() const
{
    if (!eligibleDonorCellsPtr_)
    {
        calcEligibleDonorCells();
    }

    return *eligibleDonorCellsPtr_;
}


bool Foam::oversetRegion::holePatchesPresent() const
{
    return !holeTriMesh().empty();
}


const Foam::triSurface& Foam::oversetRegion::holeTriMesh() const
{
    if (!holeTriMeshPtr_)
    {
        calcHoleTriMesh();
    }

    return *holeTriMeshPtr_;
}


const Foam::triSurfaceSearch& Foam::oversetRegion::holeSearch() const
{
    if (!holeSearchPtr_)
    {
        holeSearchPtr_ = new triSurfaceSearch
        (
            holeTriMesh()
        );
    }

    return *holeSearchPtr_;
}


const Foam::boundBox& Foam::oversetRegion::localBounds() const
{
    if (!localBoundsPtr_)
    {
        calcBounds();
    }

    return *localBoundsPtr_;
}


const Foam::boundBox& Foam::oversetRegion::globalBounds() const
{
    if (!globalBoundsPtr_)
    {
        calcBounds();
    }

    return *globalBoundsPtr_;
}


const Foam::indexedOctree<Foam::treeDataCell>&
Foam::oversetRegion::cellSearch() const
{
    if (!cellSearchPtr_)
    {
        calcCellSearch();
    }

    return *cellSearchPtr_;
}


const Foam::List<Foam::List<Foam::boundBox> >&
Foam::oversetRegion::procBoundBoxes() const
{
    if (!procBoundBoxesPtr_)
    {
        calcProcBoundBoxes();
    }

    return *procBoundBoxesPtr_;
}


void Foam::oversetRegion::update() const
{
    Info<< "oversetRegion " << name() << " update" << endl;

    fringePtr_->update();

    clearOut();
}


// ************************************************************************* //
