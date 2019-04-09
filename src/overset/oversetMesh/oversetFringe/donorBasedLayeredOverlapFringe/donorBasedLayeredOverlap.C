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
#include "faceCellsFringe.H"
#include "oversetRegion.H"
#include "addToRunTimeSelectionTable.H"

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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::donorBasedLayeredOverlapFringe::calcAddressing() const
{
    if (acceptorsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::donorBasedLayeredOverlapFringe::calcAddressing() const"
        )   << "Addressing already calculated"
            << abort(FatalError);
    }

    // Get reference to region cell zone
    const cellZone& rcz = region().zone();

    // Make a hash set to collect acceptor points
    // Note: 2 patches connecting at the corner may create a duplicate,
    // which is filtered on insertion
    labelHashSet acceptorSet;

    // Find patches and mark cells
    forAll (patchNames_, nameI)
    {
        const polyPatchID curFringePatch
        (
            patchNames_[nameI],
            mesh().boundaryMesh()
        );

        if (!curFringePatch.active())
        {
            FatalErrorIn
            (
                "void donorBasedLayeredOverlapFringe::calcAddressing() const"
            )   << "Fringe patch " << patchNames_[nameI]
                << " cannot be found"
                << abort(FatalError);
        }

        const unallocLabelList& curFaceCells =
            mesh().boundaryMesh()[curFringePatch.index()].donorBasedLayeredOverlap();

        forAll (curFaceCells, fcI)
        {
            // Check if cell is in region zone
            if (rcz.whichCell(curFaceCells[fcI]) > -1)
            {
                // Found acceptor
                acceptorSet.insert(curFaceCells[fcI]);
            }
        }
    }

    // Collect acceptors
    acceptorsPtr_ = new labelList(acceptorSet.sortedToc());

    // Holes are empty for this fringe
    fringeHolesPtr_ = new labelList();
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
        dict.lookupOrDefault<pointList>("regionCentrePoints", pointList(0))
    ),
    nLayers_(readLabel(dict.lookup("nLayers"))),
    fringeHolesPtr_(nullptr),
    acceptorsPtr_(nullptr),
    finalDonorAcceptorsPtr_(nullptr)
{
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
    forAll (connectedRegion_, crI)
    {
        // Get name of this connected region
        const word& crName = connectedRegion_[crI];

        // Find this region in the list of all regions
        const label regionID = findIndex(allRegions, crName);

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
    }

    // Sanity check: number of (optionally) specified centre points must be
    // equal to the number of connected regions
    if
    (
        regionCentrePoints_.size()
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
            "donorBasedLayeredOverlapFringe::updateIteration(donorAcceptorList&"
        )   << "finalDonorAcceptorPtr_ already allocated. Something went "
            << "wrong with the iteration procedure (flag was not updated)."
            << nl << "This should not happen for donorBasedLayeredOverlapFringe."
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
                "calcAddressing() const"
            )   << "donorBasedLayeredOverlap fringe is designed to work"
                << " with faceCells fringe as a connected region fringe."
                << nl
                << "Connected overset region " << region.name()
                << " has " << fringe.type() " fringe type. "
                << nl
                << "Proceed with care!"
                << endl;
        }

        // Update flag collecting whether all connected regions found the
        // overlap
        allFringesReady &= fringe.foundSuitableOverlap();
    }

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

            // Sanity check the validity of donor

            // Mark donors on my processor from this connected region

            // Calculate the centre for this connected region if it isn't
            // specified

            // Mark all future acceptors by moving towards the centre using
            // face-neighbour walk nLayers times

            // Now that the acceptors are marked, mark all remaining holes

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
    } // else connected fringes are not ready yet, foundSuitableOverlap flag is
      // alraedy false, so there's nothing to do

    return foundSuitablaOverlap();
}


const Foam::labelList& Foam::donorBasedLayeredOverlapFringe::fringeHoles() const
{
    if (!fringeHolesPtr_)
    {
        calcAddressing();
    }

    return *fringeHolesPtr_;
}


const Foam::labelList& Foam::donorBasedLayeredOverlapFringe::candidateAcceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
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
    if (updateFringe_)
    {
        Info<< "donorBasedLayeredOverlapFringe::update() const" << endl;

        // Clear out
        clearAddressing();
    }

    // Set flag to false and clear final donor/acceptors only
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
}


// ************************************************************************* //
