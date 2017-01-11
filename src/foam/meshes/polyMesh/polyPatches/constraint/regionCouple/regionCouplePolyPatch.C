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
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "regionCouplePolyPatch.H"
#include "polyBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "demandDrivenData.H"
#include "objectRegistry.H"
#include "polyMesh.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCouplePolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, regionCouplePolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, regionCouplePolyPatch, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::regionCouplePolyPatch::active() const
{
    // Try to find face zone
    faceZoneID zone(zoneName_, boundaryMesh().mesh().faceZones());

    if (!zone.active())
    {
        return false;
    }

    // Try to find shadow region
    if
    (
        boundaryMesh().mesh().db().parent().foundObject<polyMesh>
        (
            shadowRegionName_
        )
    )
    {
        // Shadow region present
        const polyMesh& sr = boundaryMesh().mesh().db().parent().
            objectRegistry::lookupObject<polyMesh>
            (
                shadowRegionName_
            );

        polyPatchID shadowPatch(shadowPatchName_, sr.boundaryMesh());

        // If shadow patch is active, all components are ready
        return shadowPatch.active();
    }
    else
    {
        // No shadow region
        return false;
    }
}


void Foam::regionCouplePolyPatch::calcZoneAddressing() const
{
    // Calculate patch-to-zone addressing
    if (zoneAddressingPtr_)
    {
        FatalErrorIn("void regionCouplePolyPatch::calcZoneAddressing() const")
            << "Patch to zone addressing already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "regionCouplePolyPatch::calcZoneAddressing() const for patch "
            << index() << endl;
    }

    // Calculate patch-to-zone addressing
    zoneAddressingPtr_ = new labelList(size());
    labelList& zAddr = *zoneAddressingPtr_;
    const faceZone& myZone = zone();

    for (label i = 0; i < size(); i++)
    {
        zAddr[i] = myZone.whichFace(start() + i);
    }

    // Check zone addressing
    if (zAddr.size() > 0 && min(zAddr) < 0)
    {
        FatalErrorIn("void regionCouplePolyPatch::calcZoneAddressing() const")
            << "Problem with patch-to-zone addressing for patch "
            << name()
            << ": some patch faces not found in interpolation zone"
            << abort(FatalError);
    }
}


void Foam::regionCouplePolyPatch::calcRemoteZoneAddressing() const
{
    // Calculate patch-to-zone addressing
    if (remoteZoneAddressingPtr_)
    {
        FatalErrorIn
        (
            "void regionCouplePolyPatch::calcRemoteZoneAddressing() const"
        )   << "Patch to remote zone addressing already calculated"
            << abort(FatalError);
    }

    // Once zone addressing is established, visit the opposite side and find
    // out which face data is needed for interpolation
    boolList usedShadows(shadow().zone().size(), false);

    const labelList& zAddr = zoneAddressing();

    if (master())
    {
        const labelListList& addr = patchToPatch().masterAddr();

        forAll (zAddr, mfI)
        {
            const labelList& nbrs = addr[zAddr[mfI]];

            forAll (nbrs, nbrI)
            {
                usedShadows[nbrs[nbrI]] = true;
            }
        }
    }
    else
    {
        const labelListList& addr = patchToPatch().slaveAddr();

        forAll (zAddr, mfI)
        {
            const labelList& nbrs = addr[zAddr[mfI]];

            forAll (nbrs, nbrI)
            {
                usedShadows[nbrs[nbrI]] = true;
            }
        }
    }

    // Count and pick up shadow indices
    label nShadows = 0;

    forAll (usedShadows, sI)
    {
        if (usedShadows[sI])
        {
            nShadows++;
        }
    }

    remoteZoneAddressingPtr_ = new labelList(nShadows);
    labelList& rza = *remoteZoneAddressingPtr_;

    // Reset counter for re-use
    nShadows = 0;

    forAll (usedShadows, sI)
    {
        if (usedShadows[sI])
        {
            rza[nShadows] = sI;
            nShadows++;
        }
    }
}


void Foam::regionCouplePolyPatch::calcPatchToPatch() const
{
    // Create patch-to-patch interpolation
    if (patchToPatchPtr_)
    {
        FatalErrorIn("void regionCouplePolyPatch::calcPatchToPatch() const")
            << "Patch to patch interpolation already calculated"
            << abort(FatalError);
    }

    if (master())
    {
        // Create interpolation for zones
        patchToPatchPtr_ =
            new ggiZoneInterpolation
            (
                zone()(),           // This zone reference
                shadow().zone()(),  // Shadow zone reference
                forwardT(),
                reverseT(),
                shadow().separation(), // Slave-to-master separation. Bug fix
                true,          // Patch data is complete on all processors
                SMALL,         // Non-overlapping face tolerances
                SMALL,
                true,          // Rescale weighting factors
                reject_        // Quick rejection algorithm, default BB_OCTREE
            );

        // Abort immediately if uncovered faces are present and the option
        // bridgeOverlap is not set.
        if
        (
            (
                patchToPatch().uncoveredMasterFaces().size() > 0
            && !bridgeOverlap()
            )
         || (
                patchToPatch().uncoveredSlaveFaces().size() > 0
            && !shadow().bridgeOverlap()
            )
        )
        {
            FatalErrorIn
            (
                "void regionCouplePolyPatch::calcPatchToPatch() const"
            )   << "Found uncovered faces for GGI interface "
                << name() << "/" << shadowPatchName()
                << " while the bridgeOverlap option is not set "
                << "in the boundary file." << endl
                << "This is an unrecoverable error. Aborting."
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn("void regionCouplePolyPatch::calcPatchToPatch() const")
            << "Attempting to create GGIInterpolation on a shadow"
            << abort(FatalError);
    }
}


void Foam::regionCouplePolyPatch::calcReconFaceCellCentres() const
{
    if (reconFaceCellCentresPtr_)
    {
        FatalErrorIn
        (
            "void regionCouplePolyPatch::calcReconFaceCellCentres() const"
        )   << "Reconstructed cell centres already calculated"
            << abort(FatalError);
    }

    // Create neighbouring face centres using interpolation
    if (master())
    {
        const label shadowID = shadowIndex();

        const polyMesh& sr = shadowRegion();

        // Get the transformed and interpolated shadow face cell centers
        reconFaceCellCentresPtr_ =
            new vectorField
            (
                interpolate
                (
                    sr.boundaryMesh()[shadowID].faceCellCentres()
                  - sr.boundaryMesh()[shadowID].faceCentres()
                )
              + faceCentres()
            );
    }
    else
    {
        FatalErrorIn
        (
            "void regionCouplePolyPatch::calcReconFaceCellCentres() const"
        )   << "Attempting to create reconFaceCellCentres on a shadow"
            << abort(FatalError);
    }
}


void Foam::regionCouplePolyPatch::calcLocalParallel() const
{
    // Calculate patch-to-zone addressing
    if (localParallelPtr_)
    {
        FatalErrorIn("void regionCouplePolyPatch::calcLocalParallel() const")
            << "Local parallel switch already calculated"
            << abort(FatalError);
    }

    localParallelPtr_ = new bool(false);
    bool& emptyOrComplete = *localParallelPtr_;

    // If running in parallel, all GGIs are expanded to zone size.
    // This happens on decomposition and reconstruction where
    // size and shadow size may be zero, but zone size may not
    // HJ, 1/Jun/2011
    if (!Pstream::parRun())
    {
        emptyOrComplete = true;
    }
    else
    {
        // Check that patch size is greater than the zone size.
        // This is an indication of the error where the face zone is not global
        // in a parallel run.  HJ, 9/Nov/2014
        if (size() > zone().size())
        {
            FatalErrorIn
            (
                "void regionCouplePolyPatch::calcLocalParallel() const"
            )   << "Patch size is greater than zone size for GGI patch "
                << name() << ".  This is not allowerd: "
                << "the face zone must contain all patch faces and be "
                << "global in parallel runs"
                << abort(FatalError);
        }

        // Calculate localisation on master and shadow
        emptyOrComplete =
            (
                zone().size() == size()
             && shadow().zone().size() == shadow().size()
            )
         || (size() == 0 && shadow().size() == 0);

        reduce(emptyOrComplete, andOp<bool>());
    }

    if (debug && Pstream::parRun())
    {
        Info<< "regionCouple patch Master: " << name()
            << " Slave: " << shadowPatchName() << " is ";

        if (emptyOrComplete)
        {
            Info<< "local parallel" << endl;
        }
        else
        {
            Info<< "split between multiple processors" << endl;
        }
    }
}


void Foam::regionCouplePolyPatch::calcSendReceive() const
{
    // Note: all processors will execute calcSendReceive but only master will
    // hold the information.  Therefore, pointers on slave processors
    // will remain meaningless, but for purposes of consistency
    // (of the calc-call) they will be set to zero-sized array
    // HJ, 4/Jun/2011

    if (mapPtr_)
    {
        FatalErrorIn("void regionCouplePolyPatch::calcSendReceive() const")
            << "Send-receive addressing already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "regionCouplePolyPatch::calcSendReceive() const for patch "
            << index() << endl;
    }

    if (!Pstream::parRun())
    {
        FatalErrorIn("void regionCouplePolyPatch::calcSendReceive() const")
            << "Requested calculation of send-receive addressing for a "
            << "serial run.  This is not allowed"
            << abort(FatalError);
    }

    // Gather send and receive addressing (to master)

    // Get patch-to-zone addressing
    const labelList& za = zoneAddressing();

    // Make a zone-sized field and fill it in with proc markings for processor
    // that holds and requires the data
    labelField zoneProcID(zone().size(), -1);

    forAll (za, zaI)
    {
        zoneProcID[za[zaI]] = Pstream::myProcNo();
    }

    reduce(zoneProcID, maxOp<labelField>());

    const labelList& shadowRza = shadow().remoteZoneAddressing();

    // Find out where my zone data is coming from
    labelList nRecv(Pstream::nProcs(), 0);

    // Note: only visit the data from the local zone
    forAll (shadowRza, shadowRzaI)
    {
        nRecv[zoneProcID[shadowRza[shadowRzaI]]]++;
    }

    // Make a receiving sub-map
    // It tells me which data I will receive from which processor and
    // where I need to put it into the remoteZone data before the mapping
    labelListList constructMap(Pstream::nProcs());

    // Size the receiving list
    forAll (nRecv, procI)
    {
        constructMap[procI].setSize(nRecv[procI]);
    }

    // Reset counters for processors
    nRecv = 0;

    forAll (shadowRza, shadowRzaI)
    {
        label recvProc = zoneProcID[shadowRza[shadowRzaI]];

        constructMap[recvProc][nRecv[recvProc]] = shadowRza[shadowRzaI];

        nRecv[recvProc]++;
    }

    // Make the sending sub-map
    // It tells me which data is required from me to be sent to which
    // processor

    // Algorithm
    // - expand the local zone faces with indices into a size of local zone
    // - go through remote zone addressing on all processors
    // - find out who hits my faces
    labelList localZoneIndices(zone().size(), -1);

    forAll (za, zaI)
    {
        localZoneIndices[za[zaI]] = zaI;
    }

    labelListList shadowToReceiveAddr(Pstream::nProcs());

    // Get the list of what my shadow needs to receive from my zone
    // on all other processors
    shadowToReceiveAddr[Pstream::myProcNo()] = shadowRza;
    Pstream::gatherList(shadowToReceiveAddr);
    Pstream::scatterList(shadowToReceiveAddr);

    // Now local zone indices contain the index of a local face that will
    // provide the data.  For faces that are not local, the index will be -1

    // Find out where my zone data is going to

    // Make a sending sub-map
    // It tells me which data I will send to which processor
    labelListList sendMap(Pstream::nProcs());

    // Collect local labels to be sent to each processor
    forAll (shadowToReceiveAddr, procI)
    {
        const labelList& curProcSend = shadowToReceiveAddr[procI];

        // Find out how much of my data is going to this processor
        label nProcSend = 0;

        forAll (curProcSend, sendI)
        {
            if (localZoneIndices[curProcSend[sendI]] > -1)
            {
                nProcSend++;
            }
        }

        if (nProcSend > 0)
        {
            // Collect the indices
            labelList& curSendMap = sendMap[procI];

            curSendMap.setSize(nProcSend);

            // Reset counter
            nProcSend = 0;

            forAll (curProcSend, sendI)
            {
                if (localZoneIndices[curProcSend[sendI]] > -1)
                {
                    curSendMap[nProcSend] =
                        localZoneIndices[curProcSend[sendI]];
                    nProcSend++;
                }
            }
        }
    }

    // Map will return the object of the size of remote zone
    // HJ, 9/May/2016
    mapPtr_ = new mapDistribute(zone().size(), sendMap, constructMap);
}


Foam::regionCouplePolyPatch& Foam::regionCouplePolyPatch::shadow()
{
    return const_cast<regionCouplePolyPatch&>
    (
        refCast<const regionCouplePolyPatch>
        (
            shadowRegion().boundaryMesh()[shadowIndex()]
        )
    );
}


void Foam::regionCouplePolyPatch::clearDeltas() const
{
    deleteDemandDrivenData(reconFaceCellCentresPtr_);
}


void Foam::regionCouplePolyPatch::clearGeom() const
{
    clearDeltas();

    deleteDemandDrivenData(reconFaceCellCentresPtr_);

    // Remote addressing and send-receive maps depend on the local
    // position.  Therefore, it needs to be recalculated at mesh motion.
    // Local zone addressing does not change with mesh motion
    // HJ, 23/Jun/2011
    deleteDemandDrivenData(remoteZoneAddressingPtr_);

    // localParallel depends on geometry - must be cleared!
    // HR, 11/Jul/2013
    deleteDemandDrivenData(localParallelPtr_);

    deleteDemandDrivenData(mapPtr_);
}


void Foam::regionCouplePolyPatch::clearOut() const
{
    clearGeom();

    shadowIndex_ = -1;
    zoneIndex_ = -1;

    deleteDemandDrivenData(zoneAddressingPtr_);
    deleteDemandDrivenData(patchToPatchPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCouplePolyPatch::regionCouplePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, size, start, index, bm),
    shadowRegionName_(word::null),
    shadowPatchName_(word::null),
    zoneName_(word::null),
    attached_(false),
    master_(false),
    isWall_(false),
    bridgeOverlap_(false),
    reject_(ggiZoneInterpolation::BB_OCTREE),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(NULL),
    zoneAddressingPtr_(NULL),
    remoteZoneAddressingPtr_(NULL),
    reconFaceCellCentresPtr_(NULL),
    localParallelPtr_(NULL),
    mapPtr_(NULL)
{}


Foam::regionCouplePolyPatch::regionCouplePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& shadowRegionName,
    const word& shadowPatchName,
    const word& zoneName,
    const bool attached,
    const bool master,
    const bool isWall,
    const bool bridgeOverlap,
    const ggiZoneInterpolation::quickReject reject
)
:
    coupledPolyPatch(name, size, start, index, bm),
    shadowRegionName_(shadowRegionName),
    shadowPatchName_(shadowPatchName),
    zoneName_(zoneName),
    attached_(attached),
    master_(master),
    isWall_(isWall),
    bridgeOverlap_(bridgeOverlap),
    reject_(reject),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(NULL),
    zoneAddressingPtr_(NULL),
    remoteZoneAddressingPtr_(NULL),
    reconFaceCellCentresPtr_(NULL),
    localParallelPtr_(NULL),
    mapPtr_(NULL)
{}


Foam::regionCouplePolyPatch::regionCouplePolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm),
    shadowRegionName_(dict.lookup("shadowRegion")),
    shadowPatchName_(dict.lookup("shadowPatch")),
    zoneName_(dict.lookup("zone")),
    attached_(dict.lookup("attached")),
    master_(dict.lookup("master")),
    isWall_(dict.lookup("isWall")),
    bridgeOverlap_(dict.lookup("bridgeOverlap")),
    reject_(ggiZoneInterpolation::BB_OCTREE),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(NULL),
    zoneAddressingPtr_(NULL),
    remoteZoneAddressingPtr_(NULL),
    reconFaceCellCentresPtr_(NULL),
    localParallelPtr_(NULL),
    mapPtr_(NULL)
{
    if (dict.found("quickReject"))
    {
        reject_ = ggiZoneInterpolation::quickRejectNames_.read
        (
            dict.lookup("quickReject")
        );
    }
}


Foam::regionCouplePolyPatch::regionCouplePolyPatch
(
    const regionCouplePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    shadowRegionName_(pp.shadowRegionName_),
    shadowPatchName_(pp.shadowPatchName_),
    zoneName_(pp.zoneName_),
    attached_(pp.attached_),
    master_(pp.master_),
    isWall_(pp.isWall_),
    bridgeOverlap_(pp.bridgeOverlap_),
    reject_(pp.reject_),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(NULL),
    zoneAddressingPtr_(NULL),
    remoteZoneAddressingPtr_(NULL),
    reconFaceCellCentresPtr_(NULL),
    localParallelPtr_(NULL),
    mapPtr_(NULL)
{}


Foam::regionCouplePolyPatch::regionCouplePolyPatch
(
    const regionCouplePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    shadowRegionName_(pp.shadowRegionName_),
    shadowPatchName_(pp.shadowPatchName_),
    zoneName_(pp.zoneName_),
    attached_(pp.attached_),
    master_(pp.master_),
    isWall_(pp.isWall_),
    bridgeOverlap_(pp.bridgeOverlap_),
    reject_(pp.reject_),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(NULL),
    zoneAddressingPtr_(NULL),
    remoteZoneAddressingPtr_(NULL),
    reconFaceCellCentresPtr_(NULL),
    localParallelPtr_(NULL),
    mapPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionCouplePolyPatch::~regionCouplePolyPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::regionCouplePolyPatch::coupled() const
{
    return attached_;
}


const Foam::polyMesh& Foam::regionCouplePolyPatch::shadowRegion() const
{
    if (shadowRegionName_ != Foam::word::null)
    {
        return boundaryMesh().mesh().db().parent().
            objectRegistry::lookupObject<polyMesh>
            (
                shadowRegionName()
            );
    }
    else
    {
        FatalErrorIn
        (
            "const polyMesh& regionCouplePolyPatch::shadowRegion() const"
        )   << "Requested shadowRegion which is not available"
            << abort(FatalError);

        // Dummy return
        return boundaryMesh().mesh();
    }
}


Foam::label Foam::regionCouplePolyPatch::shadowIndex() const
{
    if
    (
        shadowIndex_ == -1
     && shadowRegionName_ != Foam::word::null
     && shadowPatchName_ != Foam::word::null
    )
    {
        // Grab shadow patch index from shadow region
        const polyMesh& sr = shadowRegion();

        polyPatchID shadow(shadowPatchName_, sr.boundaryMesh());

        if (!shadow.active())
        {
            FatalErrorIn("label regionCouplePolyPatch::shadowIndex() const")
                << "Shadow patch name " << shadowPatchName_
                << " not found.  Please check your region couple "
                << "interface definition."
                << abort(FatalError);
        }

        shadowIndex_ = shadow.index();

        // Check the other side is a region couple
        if (!isA<regionCouplePolyPatch>(sr.boundaryMesh()[shadowIndex_]))
        {
            FatalErrorIn("label regionCouplePolyPatch::shadowIndex() const")
                << "Shadow of region couple patch " << name()
                << " named " << shadowPatchName()
                << " on region " << shadowRegionName()
                << " is not a region couple.  Type: "
                << boundaryMesh()[shadowIndex_].type() << nl
                << "This is not allowed.  Please check your mesh definition."
                << abort(FatalError);
        }

        // Check for region couple onto self
        if (index() == shadowIndex_ && &sr == &boundaryMesh().mesh())
        {
            FatalErrorIn("label regionCouplePolyPatch::shadowIndex() const")
                << "region couple patch " << name()
                << " created as its own shadow"
                << abort(FatalError);
        }

        // Check definition of master and slave side
        const regionCouplePolyPatch& sp =
            refCast<const regionCouplePolyPatch>
            (
                sr.boundaryMesh()[shadowIndex_]
            );

        if (master() == sp.master())
        {
            FatalErrorIn("label regionCouplePolyPatch::shadowIndex() const")
                << "Region couple patch " << name()
                << " and its shadow " << shadowPatchName()
                << " on region " << shadowRegionName()
                << ".  Clash on master-slave definition." << nl
                << "This is not allowed.  Please check your mesh definition."
                << abort(FatalError);
        }
    }

    return shadowIndex_;
}


Foam::label Foam::regionCouplePolyPatch::zoneIndex() const
{
    if (zoneIndex_ == -1 && zoneName_ != Foam::word::null)
    {
        // Grab zone patch index
        faceZoneID zone(zoneName_, boundaryMesh().mesh().faceZones());

        if (!zone.active())
        {
            FatalErrorIn("label regionCouplePolyPatch::zoneIndex() const")
                << "Face zone name " << zoneName_
                << " for region couple patch " << name()
                << " not found.  Please check your region couple "
                << "interface definition."
                << abort(FatalError);
        }

        zoneIndex_ = zone.index();
    }

    return zoneIndex_;
}


const Foam::regionCouplePolyPatch& Foam::regionCouplePolyPatch::shadow() const
{
    return refCast<const regionCouplePolyPatch>
    (
        shadowRegion().boundaryMesh()[shadowIndex()]
    );
}


const Foam::faceZone& Foam::regionCouplePolyPatch::zone() const
{
    return boundaryMesh().mesh().faceZones()[zoneIndex()];
}


void Foam::regionCouplePolyPatch::attach() const
{
    if (!attached_)
    {
        attached_ = true;
        shadow().attach();

        // Patch-to-patch interpolation does not need to be cleared,
        // only face/cell centres and interpolation factors
        // HJ, 6/Jun/2011
        //clearGeom()

        // Clear delta coefficients, but keep the rest.
        // HR, 10/Jul/2013
        clearDeltas();
    }
}


void Foam::regionCouplePolyPatch::detach() const
{
    if (attached_)
    {
        attached_ = false;
        shadow().detach();

        // Patch-to-patch interpolation does not need to be cleared,
        // only face/cell centres and interpolation factors
        // HJ, 6/Jun/2011
        //clearGeom()

        // Clear delta coefficients, but keep the rest.
        // HR, 10/Jul/2013
        clearDeltas();
    }
}


Foam::label Foam::regionCouplePolyPatch::comm() const
{
    return boundaryMesh().mesh().comm();
}


int Foam::regionCouplePolyPatch::tag() const
{
    return Pstream::msgType();
}


const Foam::labelList& Foam::regionCouplePolyPatch::zoneAddressing() const
{
    if (!zoneAddressingPtr_)
    {
        calcZoneAddressing();
    }

    return *zoneAddressingPtr_;
}


const Foam::labelList&
Foam::regionCouplePolyPatch::remoteZoneAddressing() const
{
    if (!remoteZoneAddressingPtr_)
    {
        calcRemoteZoneAddressing();
    }

    return *remoteZoneAddressingPtr_;
}


bool Foam::regionCouplePolyPatch::localParallel() const
{
    // Calculate patch-to-zone addressing
    if (!localParallelPtr_)
    {
        calcLocalParallel();
    }

    return *localParallelPtr_;
}


const Foam::ggiZoneInterpolation&
Foam::regionCouplePolyPatch::patchToPatch() const
{
    if (master())
    {
        if (!patchToPatchPtr_)
        {
            Info<< "Initializing the region couple interpolator between "
                << "master/shadow patches: "
                << name() << " and " << shadowPatchName()
                << " on region " << shadowRegionName()
                << endl;

            calcPatchToPatch();
        }

        return *patchToPatchPtr_;
    }
    else
    {
        return shadow().patchToPatch();
    }
}


const Foam::mapDistribute& Foam::regionCouplePolyPatch::map() const
{
    if (!mapPtr_)
    {
        calcSendReceive();
    }

    return *mapPtr_;
}


const Foam::vectorField&
Foam::regionCouplePolyPatch::reconFaceCellCentres() const
{
    if (!attached_)
    {
        FatalErrorIn
        (
            "const vectorField& "
            "regionCouplePolyPatch::reconFaceCellCentres() const"
        )   << "Requesting reconFaceCellCentres in detached state"
            << abort(FatalError);
    }

    if (!reconFaceCellCentresPtr_)
    {
        calcReconFaceCellCentres();
    }

    return *reconFaceCellCentresPtr_;
}


void Foam::regionCouplePolyPatch::initAddressing()
{
    if (active())
    {
        // Calculate transforms for correct GGI cut
        calcTransforms();

        if (master())
        {
            shadow().calcTransforms();
        }

        // Force zone addressing and remote zone addressing
        // (uses GGI interpolator)
        zoneAddressing();
        remoteZoneAddressing();

        // Force local parallel
        if (Pstream::parRun() && !localParallel())
        {
            // Calculate send addressing
            map();
        }
    }

    polyPatch::initAddressing();
}


void Foam::regionCouplePolyPatch::calcAddressing()
{
    polyPatch::calcAddressing();
}


void Foam::regionCouplePolyPatch::initGeometry()
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011
    if (active())
    {
        // Note: Only master calculates recon; slave uses master interpolation
        if (master())
        {
            reconFaceCellCentres();
        }
    }

    polyPatch::initGeometry();
}


void Foam::regionCouplePolyPatch::calcGeometry()
{
    polyPatch::calcGeometry();

    // Note: Calculation of transforms must be forced before the
    // reconFaceCellCentres in order to correctly set the transformation
    // in the interpolation routines.
    // HJ, 3/Jul/2009
}


void Foam::regionCouplePolyPatch::initMovePoints(const pointField& p)
{
    clearGeom();

    // Calculate transforms on mesh motion?
    calcTransforms();

    // Update interpolation for new relative position of GGI interfaces
    if (patchToPatchPtr_)
    {
        patchToPatchPtr_->movePoints
        (
            forwardT(),
            reverseT(),
            shadow().separation()
        );
    }

    // Recalculate send and receive maps
    if (active())
    {
        // Force zone addressing first
        zoneAddressing();
        remoteZoneAddressing();

        if (Pstream::parRun() && !localParallel())
        {
            // Calculate send addressing
            map();
        }
    }

    if (active() && master())
    {
        reconFaceCellCentres();
    }

    polyPatch::initMovePoints(p);
}


void Foam::regionCouplePolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);
}


void Foam::regionCouplePolyPatch::initUpdateMesh()
{
    polyPatch::initUpdateMesh();
}


void Foam::regionCouplePolyPatch::updateMesh()
{
    polyPatch::updateMesh();
    clearOut();
}


void Foam::regionCouplePolyPatch::calcTransforms() const
{
    // No transform or separation
    forwardT_.setSize(0);
    reverseT_.setSize(0);
    separation_.setSize(0);
}


void Foam::regionCouplePolyPatch::initOrder(const primitivePatch&) const
{}


bool Foam::regionCouplePolyPatch::order
(
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    // Nothing changes
    return false;
}


void Foam::regionCouplePolyPatch::syncOrder() const
{}


// Write
void Foam::regionCouplePolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("shadowRegion")
        << shadowRegionName_ << token::END_STATEMENT << nl;
    os.writeKeyword("shadowPatch")
        << shadowPatchName_ << token::END_STATEMENT << nl;
    os.writeKeyword("zone") << zoneName_
        << token::END_STATEMENT << nl;

    os.writeKeyword("attached")
        << attached_ << token::END_STATEMENT << nl;
    os.writeKeyword("master")
        << master_ << token::END_STATEMENT << nl;
    os.writeKeyword("isWall")
        << isWall_ << token::END_STATEMENT << nl;
    os.writeKeyword("bridgeOverlap") << bridgeOverlap_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
