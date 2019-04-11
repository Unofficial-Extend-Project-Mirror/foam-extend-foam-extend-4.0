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

#include "donorBasedLayeredOverlapFringe.H"
#include "oversetMesh.H"
#include "oversetRegion.H"
#include "faceCellsFringe.H"
#include "oversetRegion.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(donorBasedLayeredOverlapFringe, 0);
    addToRunTimeSelectionTable
    (
        oversetFringe,
        donorBasedLayeredOverlapFringe,
        dictionary
    );
}


const Foam::debug::tolerancesSwitch
Foam::donorBasedLayeredOverlapFringe::distTol_
(
    "donorBasedLayeredOverlapDistanceTolerance",
    0.0
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::donorBasedLayeredOverlapFringe::calcAddressing() const
{
    // Make sure that either acceptorsPtr is unnalocated or if it is allocated,
    // that it is empty
    if (acceptorsPtr_ && !acceptorsPtr_->empty())
    {
        FatalErrorIn
        (
            "void Foam::donorBasedLayeredOverlapFringe::calcAddressing() const"
        )   << "Addressing already calculated"
            << abort(FatalError);
    }

    // Get list of all overset regions
    const PtrList<oversetRegion>& allRegions =
        this->region().overset().regions();

    // Loop through all connected regions and check whether the fringe overlap
    // has been found for all of them
    bool allFringesReady = true;
    forAll (connectedRegionIDs_, crI)
    {
        // Get ID of this region
        const label& regionID = connectedRegionIDs_[crI];

        // Get the overset region
        const oversetRegion& region = allRegions[regionID];

        // Get the fringe algorithm from the region
        const oversetFringe& fringe = region.fringe();

        // If this is not faceCells fringe, issue a Warning. This fringe
        // selection algorithm is intended to work only with faceCells fringe on
        // the other side. VV, 9/Apr/2019
        if (!isA<faceCellsFringe>(fringe))
        {
            WarningIn
            (
                "void Foam::donorBasedLayeredOverlapFringe::"
                "updateIteration(donorAcceptorList&) const"
            )   << "donorBasedLayeredOverlap fringe is designed to work"
                << " with faceCells fringe as a connected region fringe."
                << nl
                << "Connected overset region " << region.name()
                << " has " << fringe.type() << " fringe type. "
                << nl
                << "Proceed with care!"
                << endl;
        }

        // Update flag collecting whether all connected regions found the
        // overlap
        allFringesReady &= fringe.foundSuitableOverlap();
    }

    // Sets containing all acceptors and all holes for all connected regions
    const polyMesh& mesh = this->mesh();
    labelHashSet allAcceptors(0.02*mesh.nCells());
    labelHashSet allFringeHoles(0.02*mesh.nCells());

    if (allFringesReady)
    {
        // Loop through connected regions
        forAll (connectedRegionIDs_, crI)
        {
            // Get ID of this region
            const label& regionID = connectedRegionIDs_[crI];

            // Get fringe of the connected region
            const oversetFringe& fringe = allRegions[regionID].fringe();

            // The fringe should be finalized, which means we may take a const
            // reference to its final donor acceptors
            const donorAcceptorList& crDonorAcceptorPairs =
                fringe.finalDonorAcceptors();

            // Hash set containing donors
            labelHashSet donors(6*crDonorAcceptorPairs.size());

            // Initialize centre of the donors of this connected region in order
            // to search in a given direction
            vector centrePoint(vector::zero);

            // Loop through all donor/acceptors
            forAll (crDonorAcceptorPairs, daI)
            {
                // Get this donor/acceptor pair
                const donorAcceptor& daPair = crDonorAcceptorPairs[daI];

                // Check whether all donors have been found
                if (!daPair.donorFound())
                {
                    FatalErrorIn
                    (
                        "donorBasedLayeredOverlapFringe::"
                        "updateIteration(donorAcceptorList&) const"
                    )   << "Donor not found for donor/acceptor pair " << daI
                        << nl
                        << "Donor/acceptor data: " << daPair
                        << nl
                        << "In connected region: " << allRegions[regionID].name()
                        << abort(FatalError);
                }

                // Mark donors on my processor from this connected region. Note
                // that the check has been made in constructor to make sure that
                // this region is the only donor region for the connected region
                if (daPair.donorProcNo() == Pstream::myProcNo())
                {
                    // Get donor index
                    const label& dI = daPair.donorCell();

                    // Insert donor into the hash set
                    if (donors.insert(dI))
                    {
                        // Donor has been inserted (not previously found in the
                        // hash set), add donor point to centre point (the
                        // centre point will be calculated later on as
                        // arithmetic mean)
                        centrePoint += daPair.donorPoint();
                    }

                    // Loop through extended donor cells
                    const donorAcceptor::DynamicLabelList& extDonors =
                        daPair.extendedDonorCells();
                    const donorAcceptor::DynamicPointList& extDonorPoints =
                        daPair.extendedDonorPoints();

                    forAll (extDonors, i)
                    {
                        // Get donor index
                        const label& edI = extDonors[i];

                        // Inser extended donor into the hash set
                        if (donors.insert(edI))
                        {
                            // Donor has been inserted (not previously found in
                            // the hash set), add extended donor point as well
                            centrePoint += extDonorPoints[i];
                        }
                    } // End for all extended donors
                } // End if this donor is on my processor
            } // End for all (master) donor cells

            // Use the centre point as specified by the user if it was specified
            // (if the regionCentrePoints_ list is not empty). This avoids
            // parallel communication as well.
            if (!regionCentrePoints_.empty())
            {
                // Use specified centre point, discarding the data we calculated
                // above
                centrePoint = regionCentrePoints_[crI];
            }
            else
            {
                // User did not specify centre points and the centre point holds
                // the sum of all the points. Reduce the data
                reduce(centrePoint, sumOp<vector>());
                const label nUniqueDonors =
                    returnReduce(donors.size(), sumOp<label>());

                // Calculate the final centre point by finding the arithmetic mean
                centrePoint /= nUniqueDonors;
            }

            // We now have a collection of all donors for this connected region
            // and the centre point to move to. Let's collect the acceptors
            labelHashSet acceptors(donors.size()); // Reasonable size estimate

            // Get necessary mesh data (from polyMesh/primitiveMesh)
            const vectorField& cc = mesh.cellCentres();
            const vectorField& fc = mesh.faceCentres();
            const cellList& meshCells = mesh.cells();
            const unallocLabelList& owner = mesh.faceOwner();
            const unallocLabelList& neighbour = mesh.faceNeighbour();

            // Mark cells that are eligible to be acceptors (not donors)
            boolList eligibleCells(mesh.nCells(), true);
            forAllConstIter (labelHashSet, donors, iter)
            {
                eligibleCells[iter.key()] = false;
            }

            // Loop nLayers away from initial donors
            for (label i = 0; i < nLayers_; ++i)
            {
                // Face markup for propagation
                boolList propagateFace(mesh.nFaces(), false);

                // Loop through all donors and mark faces that are pointing
                // towards the centre point and have an eligible neighbour
                forAllConstIter (labelHashSet, donors, iter)
                {
                    // Get the cell index and the cell
                    const label& cellI = iter.key();
                    const cell& cFaces = meshCells[cellI];

                    // Get cell centre of this donor and calculate distance to
                    // centre point
                    const vector& donorCentre = cc[cellI];
                    const scalar donorCentreToRegionCentreDist =
                        mag(donorCentre - centrePoint);

                    // Loop through all faces of the cell
                    forAll (cFaces, i)
                    {
                        // Get face index (global)
                        const label& faceI = cFaces[i];

                        // Get face centre and calculate distance to centre
                        // point
                        const vector& faceCentre = fc[faceI];
                        const scalar faceCentreToRegionCentreDist =
                            mag(faceCentre - centrePoint);

                        if
                        (
                            faceCentreToRegionCentreDist
                          - donorCentreToRegionCentreDist
                          < distTol_
                        )
                        {
                            // Face is closer to the centre point than cell: we
                            // are moving in the right direction. Mark the face
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

                // Loop through all faces and append acceptors
                for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
                {
                    if (propagateFace[faceI])
                    {
                        // Face is marked, select owner or neighbour
                        const label& own = owner[faceI];
                        const label& nei = neighbour[faceI];

                        // Either owner or neighbour may be eligible, not both
                        if (eligibleCells[own])
                        {
                            // Owner cell is not a donor, insert it
                            acceptors.insert(own);

                            // Mark as ineligible
                            eligibleCells[own] = false;
                        }
                        else if (eligibleCells[nei])
                        {
                            // Neighbour cell is not a donor, insert it
                            acceptors.insert(nei);

                            // Mark as ineligible
                            eligibleCells[nei] = false;
                        }
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

                        if (eligibleCells[own])
                        {
                            // Face cell is not a donor, insert it
                            acceptors.insert(own);

                            // Mark as ineligible
                            eligibleCells[own] = false;
                        }
                    }
                }

                // Special treatment for last iteration
                if (i == nLayers_ - 1)
                {
                    // This is not the last iteration, transfer acceptors into
                    // donors
                    donors.transfer(acceptors);

                    // Resize acceptors list
                    acceptors.resize(donors.size());
                }
            } // End for specified number of layers

            // At this point, we have the final set of acceptors and we marked
            // all cells that are ineligible (either donor or acceptor). The
            // remaining thing to do is to mark the interior holes

            // Create a hash set containing fringe holes, initialized with
            // acceptors in order to avoid having special conditions in the loop
            // below for the first pass
            labelHashSet fringeHoles(10*acceptors.size());

            // Collect holes until there are no holes to collect
            label nAddedHoles;
            do
            {
                // Reset number of newly added holes
                nAddedHoles = 0;

                // Face markup for propagation
                boolList propagateFace(mesh.nFaces(), false);

                // Loop through all acceptors and mark faces that are pointing
                // towards the centre point and have eligible neighbour (not an
                // acceptor or donor)
                forAllConstIter (labelHashSet, fringeHoles, iter)
                {
                    // Get the cell index and the cell
                    const label& cellI = iter.key();
                    const cell& cFaces = meshCells[cellI];

                    // Note: there's no need to check for the distance here
                    // because there's always at least one "buffer" layer
                    // towards the outer side that consists of donors, which are
                    // marked as ineligible at the beginning

                    // Loop through all faces of the cell and mark all of them
                    // for propagation
                    forAll (cFaces, i)
                    {
                        propagateFace[cFaces[i]] = true;
                    }
                }

                // Sync the face list across processor boundaries
                syncTools::syncFaceList
                (
                    mesh,
                    propagateFace,
                    orOp<bool>(),
                    false
                );

                // Loop through all faces and append interior holes
                for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
                {
                    if (propagateFace[faceI])
                    {
                        // Face is marked, select owner or neighbour
                        const label& own = owner[faceI];
                        const label& nei = neighbour[faceI];

                        // Either owner or neighbour may be eligible, not both
                        if (eligibleCells[own])
                        {
                            // Owner cell is not a hole (or an acceptor in the
                            // first iteration), insert it
                            if (fringeHoles.insert(own))
                            {
                                // Count number of added holes in this iteration
                                ++nAddedHoles;
                            }

                            // Mark as ineligible
                            eligibleCells[own] = false;
                        }
                        else if (eligibleCells[nei])
                        {
                            // Neighbour cell is not a hole (or an acceptor in
                            // the first iteration), insert it
                            if (fringeHoles.insert(nei))
                            {
                                // Count number of added holes in this iteration
                                ++nAddedHoles;
                            }

                            // Mark as ineligible
                            eligibleCells[nei] = false;
                        }
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

                        if (eligibleCells[own])
                        {
                            // Face cell is not a hole (or an acceptor in the
                            // first iteration), inser it
                            if (fringeHoles.insert(own))
                            {
                                ++nAddedHoles;
                            }

                            // Mark as ineligible
                            eligibleCells[own] = false;
                        }
                    }
                }

                // We moved one layer "inside" the fringe. Keep going until
                // there are no more holes to add

            } while (nAddedHoles != 0);

            // Finally, we have collected all the fringe holes and acceptors
            // from this connected region. Append them to the global sets
            allAcceptors += acceptors;
            allFringeHoles += fringeHoles;

        } // End for all connected regions

        // Set acceptors and holes from the data for all regions
        acceptorsPtr_ = new labelList(allAcceptors.sortedToc());
        fringeHolesPtr_ = new labelList(allFringeHoles.sortedToc());
    }
    else
    {
        // Connected fringes are not ready, allocate empty lists for acceptors
        // and holes, which will be deleted when asked for again from the
        // iterative procedure (see candidateAcceptors() and fringeHoles())
        acceptorsPtr_ = new labelList(0);
        fringeHolesPtr_ = new labelList(0);
    }
}


void Foam::donorBasedLayeredOverlapFringe::clearAddressing() const
{
    deleteDemandDrivenData(fringeHolesPtr_);
    deleteDemandDrivenData(acceptorsPtr_);
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::donorBasedLayeredOverlapFringe::donorBasedLayeredOverlapFringe
(
    const fvMesh& mesh,
    const oversetRegion& region,
    const dictionary& dict
)
:
    oversetFringe(mesh, region, dict),
    connectedRegionIDs_(),
    regionCentrePoints_
    (
        dict.lookupOrDefault<List<point> >
        (
            "regionCentrePoints",
            List<point>(0)
        )
    ),
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
            "donorBasedLayeredOverlapFringe::"
            "donorBasedLayeredOverlapFringe\n"
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

    // Read names of connected regions
    const wordList connectedRegionNames(dict.lookup("connectedRegions"));

    // Set size of the list containing IDs
    connectedRegionIDs_.setSize(connectedRegionNames.size());

    // Get list of all overset regions
    const PtrList<oversetRegion>& allRegions =
        this->region().overset().regions();

    // Create list of all region names for easy lookup
    wordList allRegionNames(allRegions.size());
    forAll (allRegionNames, arI)
    {
        allRegionNames[arI] = allRegions[arI].name();
    }

    // Loop through all regions and check whether the overlap has been found
    forAll (connectedRegionNames, crI)
    {
        // Get name of this connected region
        const word& crName = connectedRegionNames[crI];

        // Find this region in the list of all regions
        const label regionID = findIndex(allRegionNames, crName);

        if (regionID == -1)
        {
            FatalErrorIn
            (
                "donorBasedLayeredOverlapFringe::"
                "donorBasedLayeredOverlapFringe\n"
                "(\n"
                "    const fvMesh& mesh,\n"
                "    const oversetRegion& region,\n"
                "    const dictionary& dict\n"
                ")"
            )   << "Region " << crName << " not found in list of regions."
                << "List of overset regions: " << allRegionNames
                << abort(FatalError);
        }

        // Collect the region index in the list
        connectedRegionIDs_[crI] = regionID;

        // Sanity check: if the specified connected donor region has more than 1
        // donor regions, this fringe algorithm is attempted to be used for
        // something that's not intended. Issue an error
        if (allRegions[regionID].donorRegions().size() != 1)
        {
            FatalErrorIn
            (
                "donorBasedLayeredOverlapFringe::"
                "donorBasedLayeredOverlapFringe\n"
                "(\n"
                "    const fvMesh& mesh,\n"
                "    const oversetRegion& region,\n"
                "    const dictionary& dict\n"
                ")"
            )   << "Region " << crName << " specified as connected region, but"
                << " that region has "
                << allRegions[regionID].donorRegions().size() << " donor regions."
                << abort(FatalError);
        }

        // Sanity check whether the donor region of connected region is actually
        // this region
        if (regionID != this->region().index())
        {
            FatalErrorIn
            (
                "donorBasedLayeredOverlapFringe::"
                "donorBasedLayeredOverlapFringe\n"
                "(\n"
                "    const fvMesh& mesh,\n"
                "    const oversetRegion& region,\n"
                "    const dictionary& dict\n"
                ")"
            )   << "The donor region of region " << crName
                << " should be only region " << this->region().name()
                << abort(FatalError);
        }
    }

    // Sanity check: number of (optionally) specified centre points must be
    // equal to the number of connected regions
    if
    (
        !regionCentrePoints_.empty()
     && (regionCentrePoints_.size() != connectedRegionIDs_.size())
    )
    {
        // The list is not empty and the size of the list is not the same as the
        // size of the connected regions. This is a problem
        FatalErrorIn
        (
            "donorBasedLayeredOverlapFringe::"
            "donorBasedLayeredOverlapFringe\n"
            "(\n"
            "    const fvMesh& mesh,\n"
            "    const oversetRegion& region,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "You have specified "
            << regionCentrePoints_.size()
            << " regionCentrePoints, while specifying "
            << connectedRegionIDs_.size()
            << " connectedRegions."
            << nl
            << "If you'd like to avoid using automatic centre point detection,"
            << " make sure to specify centre points for all connected regions."
            << abort(FatalError);
    }

    // Sanity check: if the user has specified centre points, dynamic mesh
    // simulations may be problematic. Issue a warning
    if
    (
        !regionCentrePoints_.empty()
     && isA<dynamicFvMesh>(this->mesh())
    )
    {
        WarningIn
        (
            "donorBasedLayeredOverlapFringe::"
            "donorBasedLayeredOverlapFringe\n"
            "(\n"
            "    const fvMesh& mesh,\n"
            "    const oversetRegion& region,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "You have specified regionCentrePoints for a dynamic mesh"
            << " simulation."
            << nl
            << "Make sure that the centre points always remain in/near the"
            << " centre of donors in connected regions!"
            << nl
            << "Proceed with care!"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::donorBasedLayeredOverlapFringe::~donorBasedLayeredOverlapFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::donorBasedLayeredOverlapFringe::updateIteration
(
    donorAcceptorList& donorAcceptorRegionData
) const
{
    // If the donorAcceptor list has been allocated, something went wrong with
    // the iteration procedure (not-updated flag): this function has been called
    // more than once, which should not happen for donorBasedLayeredOverlapFringe
    if (finalDonorAcceptorsPtr_)
    {
        FatalErrorIn
        (
            "donorBasedLayeredOverlapFringe::"
            "updateIteration(donorAcceptorList&) const"
        )   << "finalDonorAcceptorPtr_ already allocated. Something went "
            << "wrong with the iteration procedure (flag was not updated)."
            << nl << "This should not happen for donorBasedLayeredOverlapFringe."
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


const Foam::labelList& Foam::donorBasedLayeredOverlapFringe::fringeHoles() const
{
    if (!fringeHolesPtr_)
    {
        calcAddressing();
    }
    else if (fringeHolesPtr_->empty())
    {
        // Fringe holes pointer is empty, delete both acceptorPtr_ and
        // fringeHolesPtr_ and calculate addressing. Since this fringe strategy
        // depends on other fringes, we need to have a special control in the
        // iterative algorithm that will start only when all the others are done
        // See calcAddressing() for details
        deleteDemandDrivenData(fringeHolesPtr_);
        deleteDemandDrivenData(acceptorsPtr_);
    }

    return *fringeHolesPtr_;
}


const Foam::labelList& Foam::donorBasedLayeredOverlapFringe::candidateAcceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }
    else if (acceptorsPtr_->empty())
    {
        // Acceptors pointer is empty, delete both acceptorPtr_ and
        // fringeHolesPtr_ and calculate addressing. Since this fringe strategy
        // depends on other fringes, we need to have a special control in the
        // iterative algorithm that will start only when all the others are done
        // See calcAddressing() for details
        deleteDemandDrivenData(fringeHolesPtr_);
        deleteDemandDrivenData(acceptorsPtr_);
    }

    return *acceptorsPtr_;
}


Foam::donorAcceptorList& Foam::donorBasedLayeredOverlapFringe::finalDonorAcceptors() const
{
    if (!finalDonorAcceptorsPtr_)
    {
        FatalErrorIn("donorBasedLayeredOverlapFringe::finalDonorAcceptors()")
            << "finalDonorAcceptorPtr_ not allocated. Make sure you have "
            << "called donorBasedLayeredOverlapFringe::updateIteration() before asking for "
            << "final set of donor/acceptor pairs."
            << abort(FatalError);
    }

    if (!foundSuitableOverlap())
    {
        FatalErrorIn("donorBasedLayeredOverlapFringe::finalDonorAcceptors()")
            << "Attemted to access finalDonorAcceptors but suitable overlap "
            << "has not been found. This is not allowed. "
            << abort(FatalError);
    }

    return *finalDonorAcceptorsPtr_;
}


void Foam::donorBasedLayeredOverlapFringe::update() const
{
    Info<< "donorBasedLayeredOverlapFringe::update() const" << endl;

    // Clear out
    clearAddressing();

    // Set flag to false and clear final donor/acceptors only
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
    updateSuitableOverlapFlag(false);
}


// ************************************************************************* //
