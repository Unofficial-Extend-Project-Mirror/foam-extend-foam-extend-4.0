/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Contributor
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "ggiPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "demandDrivenData.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"
#include "SubField.H"
#include "foamTime.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebugWithDescription
    (
        ggiPolyPatch,
        0,
        "If value > 1, write uncovered GGI patch facets to VTK file"
    );

    addToRunTimeSelectionTable(polyPatch, ggiPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, ggiPolyPatch, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::ggiPolyPatch::active() const
{
    polyPatchID shadow(shadowName_, boundaryMesh());
    faceZoneID zone(zoneName_, boundaryMesh().mesh().faceZones());

    return shadow.active() && zone.active();
}


void Foam::ggiPolyPatch::calcZoneAddressing() const
{
    // Calculate patch-to-zone addressing
    if (zoneAddressingPtr_)
    {
        FatalErrorIn("void ggiPolyPatch::calcZoneAddressing() const")
            << "Patch to zone addressing already calculated"
            << abort(FatalError);
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
        Info<< "myZone: " << myZone << nl
            << "my start and size: " << start() << " and " << size() << nl
            << "zAddr: " << zAddr << endl;

        FatalErrorIn("void ggiPolyPatch::calcZoneAddressing() const")
            << "Problem with patch-to-zone addressing: some patch faces "
            << "not found in interpolation zone"
            << abort(FatalError);
    }
}


void Foam::ggiPolyPatch::calcRemoteZoneAddressing() const
{
    // Calculate patch-to-zone addressing
    if (remoteZoneAddressingPtr_)
    {
        FatalErrorIn("void ggiPolyPatch::calcRemoteZoneAddressing() const")
            << "Patch to remote zone addressing already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "ggiPolyPatch::calcRemoteZoneAddressing() const for patch "
            << name() << endl;
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


void Foam::ggiPolyPatch::calcPatchToPatch() const
{
    // Create patch-to-patch interpolation
    if (patchToPatchPtr_)
    {
        FatalErrorIn("void ggiPolyPatch::calcPatchToPatch() const")
            << "Patch to patch interpolation already calculated"
            << abort(FatalError);
    }

    if (master())
    {
        if (debug)
        {
            InfoIn("void ggiPolyPatch::calcPatchToPatch() const")
                << "Calculating patch to patch interpolation for patch"
                << name() << endl;
        }

        // Create interpolation for zones
        patchToPatchPtr_ =
            new ggiZoneInterpolation
            (
                zone()(),           // This zone reference
                shadow().zone()(),  // Shadow zone reference
                forwardT(),
                reverseT(),
                -separation(), // Slave-to-master separation: Use - localValue
                // Bug fix, delayed slave evaluation causes error
                // HJ, 30/Jun/2013
                0,             // Non-overlapping face tolerances
                0,             // HJ, 24/Oct/2008
                true,          // Rescale weighting factors.  Bug fix, MB.
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
            FatalErrorIn("void ggiPolyPatch::calcPatchToPatch() const")
                << "Found uncovered faces for GGI interface "
                << name() << "/" << shadowName()
                << " while the bridgeOverlap option is not set "
                << "in the boundary file." << endl
                << "This is an unrecoverable error. Aborting."
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn("void ggiPolyPatch::calcPatchToPatch() const")
            << "Attempting to create GGIInterpolation on a shadow"
            << abort(FatalError);
    }
}


void Foam::ggiPolyPatch::calcReconFaceCellCentres() const
{
    if (reconFaceCellCentresPtr_)
    {
        FatalErrorIn
        (
            "void ggiPolyPatch::calcReconFaceCellCentres() const"
        )   << "Reconstructed cell centres already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        InfoIn("void ggiPolyPatch::calcReconFaceCellCentres() const")
            << "Calculating recon centres for patch"
            << name() << endl;
    }

    // Create neighbouring face centres using interpolation
    if (master())
    {
        const label shadowID = shadowIndex();

        // Get the transformed and interpolated shadow face cell centers
        reconFaceCellCentresPtr_ =
            new vectorField
            (
                interpolate
                (
                    boundaryMesh()[shadowID].faceCellCentres()
                  - boundaryMesh()[shadowID].faceCentres()
                )
              + faceCentres()
            );
    }
    else
    {
        FatalErrorIn("void ggiPolyPatch::calcReconFaceCellCentres() const")
            << "Attempting to create reconFaceCellCentres on a shadow"
            << abort(FatalError);
    }
}


void Foam::ggiPolyPatch::calcLocalParallel() const
{
    // Calculate patch-to-zone addressing
    if (localParallelPtr_)
    {
        FatalErrorIn("void ggiPolyPatch::calcLocalParallel() const")
            << "Local parallel switch already calculated"
            << abort(FatalError);
    }

    localParallelPtr_ = new bool(false);
    bool& emptyOrComplete = *localParallelPtr_;

    // If running in serial, all GGIs are expanded to zone size.
    // This happens on decomposition and reconstruction where
    // size and shadow size may be zero, but zone size may not
    // HJ, 1/Jun/2011
    if (!Pstream::parRun())
    {
        emptyOrComplete = false;
    }
    else
    {
        // Check that patch size is greater than the zone size.
        // This is an indication of the error where the face zone is not global
        // in a parallel run.  HJ, 9/Nov/2014
        if (size() > zone().size())
        {
            FatalErrorIn("void ggiPolyPatch::calcLocalParallel() const")
                << "Patch size is greater than zone size for GGI patch "
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
        Info<< "GGI patch Master: " << name()
            << " Slave: " << shadowName() << " is ";

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


void Foam::ggiPolyPatch::calcSendReceive() const
{
    // Note: all processors will execute calcSendReceive but only master will
    // hold the information.  Therefore, pointers on slave processors
    // will remain meaningless, but for purposes of consistency
    // (of the calc-call) they will be set to zero-sized array
    // HJ, 4/Jun/2011

    if (receiveAddrPtr_ || sendAddrPtr_)
    {
        FatalErrorIn("void ggiPolyPatch::calcSendReceive() const")
            << "Send-receive addressing already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "ggiPolyPatch::calcSendReceive() const for patch "
            << name() << endl;
    }

    if (!Pstream::parRun())
    {
        FatalErrorIn("void ggiPolyPatch::calcSendReceive() const")
            << "Requested calculation of send-receive addressing for a "
            << "serial run.  This is not allowed"
            << abort(FatalError);
    }

    // Master will receive and store the maps
    if (Pstream::master())
    {
        receiveAddrPtr_ = new labelListList(Pstream::nProcs());
        labelListList& rAddr = *receiveAddrPtr_;

        sendAddrPtr_ = new labelListList(Pstream::nProcs());
        labelListList& sAddr = *sendAddrPtr_;

        // Insert master
        rAddr[0] = zoneAddressing();

        for (label procI = 1; procI < Pstream::nProcs(); procI++)
        {
            // Note: must use normal comms because the size of the
            // communicated lists is unknown on the receiving side
            // HJ, 4/Jun/2011

            // Opt: reconsider mode of communication
            IPstream ip(Pstream::scheduled, procI);

            rAddr[procI] = labelList(ip);

            sAddr[procI] = labelList(ip);
        }
    }
    else
    {
        // Create dummy pointers: only master processor stores maps
        receiveAddrPtr_ = new labelListList();
        sendAddrPtr_ = new labelListList();

        // Send information to master
        const labelList& za = zoneAddressing();
        const labelList& ra = remoteZoneAddressing();

        // Note: must use normal comms because the size of the
        // communicated lists is unknown on the receiving side
        // HJ, 4/Jun/2011

        // Opt: reconsider mode of communication
        OPstream op(Pstream::scheduled, Pstream::masterNo());

        // Send local and remote addressing to master
        op << za << ra;
    }
}


const Foam::labelListList& Foam::ggiPolyPatch::receiveAddr() const
{
    if (!receiveAddrPtr_)
    {
        calcSendReceive();
    }

    return *receiveAddrPtr_;
}


const Foam::labelListList& Foam::ggiPolyPatch::sendAddr() const
{
    if (!sendAddrPtr_)
    {
        calcSendReceive();
    }

    return *sendAddrPtr_;
}


void Foam::ggiPolyPatch::clearGeom()
{
    deleteDemandDrivenData(reconFaceCellCentresPtr_);

    // Remote addressing and send-receive maps depend on the local
    // position.  Therefore, it needs to be recalculated at mesh motion.
    // Local zone addressing does not change with mesh motion
    // HJ, 23/Jun/2011
    deleteDemandDrivenData(remoteZoneAddressingPtr_);

    deleteDemandDrivenData(receiveAddrPtr_);
    deleteDemandDrivenData(sendAddrPtr_);
}


void Foam::ggiPolyPatch::clearOut()
{
    clearGeom();

    shadowIndex_ = -1;
    zoneIndex_ = -1;

    deleteDemandDrivenData(zoneAddressingPtr_);
    deleteDemandDrivenData(patchToPatchPtr_);
    deleteDemandDrivenData(localParallelPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ggiPolyPatch::ggiPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, size, start, index, bm),
    shadowName_("initializeMe"),
    zoneName_("initializeMe"),
    bridgeOverlap_(false),
    reject_(ggiZoneInterpolation::BB_OCTREE),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(NULL),
    zoneAddressingPtr_(NULL),
    remoteZoneAddressingPtr_(NULL),
    reconFaceCellCentresPtr_(NULL),
    localParallelPtr_(NULL),
    receiveAddrPtr_(NULL),
    sendAddrPtr_(NULL)
{}


Foam::ggiPolyPatch::ggiPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& shadowName,
    const word& zoneName,
    const bool bridgeOverlap,
    const ggiZoneInterpolation::quickReject reject
)
:
    coupledPolyPatch(name, size, start, index, bm),
    shadowName_(shadowName),
    zoneName_(zoneName),
    bridgeOverlap_(bridgeOverlap),
    reject_(reject),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(NULL),
    zoneAddressingPtr_(NULL),
    remoteZoneAddressingPtr_(NULL),
    reconFaceCellCentresPtr_(NULL),
    localParallelPtr_(NULL),
    receiveAddrPtr_(NULL),
    sendAddrPtr_(NULL)
{}


Foam::ggiPolyPatch::ggiPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm),
    shadowName_(dict.lookup("shadowPatch")),
    zoneName_(dict.lookup("zone")),
    bridgeOverlap_(dict.lookup("bridgeOverlap")),
    reject_(ggiZoneInterpolation::BB_OCTREE),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(NULL),
    zoneAddressingPtr_(NULL),
    remoteZoneAddressingPtr_(NULL),
    reconFaceCellCentresPtr_(NULL),
    localParallelPtr_(NULL),
    receiveAddrPtr_(NULL),
    sendAddrPtr_(NULL)
{
    if (dict.found("quickReject"))
    {
        reject_ = ggiZoneInterpolation::quickRejectNames_.read
        (
            dict.lookup("quickReject")
        );
    }
}


Foam::ggiPolyPatch::ggiPolyPatch
(
    const ggiPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    shadowName_(pp.shadowName_),
    zoneName_(pp.zoneName_),
    bridgeOverlap_(pp.bridgeOverlap_),
    reject_(pp.reject_),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(NULL),
    zoneAddressingPtr_(NULL),
    remoteZoneAddressingPtr_(NULL),
    reconFaceCellCentresPtr_(NULL),
    localParallelPtr_(NULL),
    receiveAddrPtr_(NULL),
    sendAddrPtr_(NULL)
{}


//- Construct as copy, resetting the face list and boundary mesh data
Foam::ggiPolyPatch::ggiPolyPatch
(
    const ggiPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    shadowName_(pp.shadowName_),
    zoneName_(pp.zoneName_),
    bridgeOverlap_(pp.bridgeOverlap_),
    reject_(pp.reject_),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(NULL),
    zoneAddressingPtr_(NULL),
    remoteZoneAddressingPtr_(NULL),
    reconFaceCellCentresPtr_(NULL),
    localParallelPtr_(NULL),
    receiveAddrPtr_(NULL),
    sendAddrPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ggiPolyPatch::~ggiPolyPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::ggiPolyPatch::shadowIndex() const
{
    if (shadowIndex_ == -1 && shadowName_ != Foam::word::null)
    {
        // Grab shadow patch index
        polyPatchID shadow(shadowName_, boundaryMesh());

        if (!shadow.active())
        {
            FatalErrorIn("label ggiPolyPatch::shadowIndex() const")
                << "Shadow patch name " << shadowName_
                << " not found.  Please check your GGI interface definition."
                << abort(FatalError);
        }

        shadowIndex_ = shadow.index();

        // Check the other side is a ggi
        if (!isA<ggiPolyPatch>(boundaryMesh()[shadowIndex_]))
        {
            FatalErrorIn("label ggiPolyPatch::shadowIndex() const")
                << "Shadow of ggi patch " << name()
                << " named " << shadowName() << " is not a ggi.  Type: "
                << boundaryMesh()[shadowIndex_].type() << nl
                << "This is not allowed.  Please check your mesh definition."
                << abort(FatalError);
        }

        // Check for GGI onto self
        if (index() == shadowIndex_)
        {
            FatalErrorIn("label ggiPolyPatch::shadowIndex() const")
                << "ggi patch " << name() << " created as its own shadow"
                << abort(FatalError);
        }
    }

    return shadowIndex_;
}


Foam::label Foam::ggiPolyPatch::zoneIndex() const
{
    if (zoneIndex_ == -1 && zoneName_ != Foam::word::null)
    {
        // Grab zone patch index
        faceZoneID zone(zoneName_, boundaryMesh().mesh().faceZones());

        if (!zone.active())
        {
            FatalErrorIn("label ggiPolyPatch::zoneIndex() const")
                << "Face zone name " << zoneName_
                << " for GGI patch " << name()
                << " not found.  Please check your GGI interface definition."
                << abort(FatalError);
        }

        zoneIndex_ = zone.index();
    }

    return zoneIndex_;
}


const Foam::ggiPolyPatch& Foam::ggiPolyPatch::shadow() const
{
    return refCast<const ggiPolyPatch>(boundaryMesh()[shadowIndex()]);
}


const Foam::faceZone& Foam::ggiPolyPatch::zone() const
{
    return boundaryMesh().mesh().faceZones()[zoneIndex()];
}


const Foam::labelList& Foam::ggiPolyPatch::zoneAddressing() const
{
    if (!zoneAddressingPtr_)
    {
        calcZoneAddressing();
    }

    return *zoneAddressingPtr_;
}


const Foam::labelList& Foam::ggiPolyPatch::remoteZoneAddressing() const
{
    if (!remoteZoneAddressingPtr_)
    {
        calcRemoteZoneAddressing();
    }

    return *remoteZoneAddressingPtr_;
}


bool Foam::ggiPolyPatch::localParallel() const
{
    // Calculate patch-to-zone addressing
    if (!localParallelPtr_)
    {
        calcLocalParallel();
    }

    return *localParallelPtr_;
}


const Foam::ggiZoneInterpolation& Foam::ggiPolyPatch::patchToPatch() const
{
    if (master())
    {
        if (!patchToPatchPtr_)
        {
            Info<< "Initializing the GGI interpolator between "
                << "master/shadow patches: "
                << name() << "/" << shadowName()
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


const Foam::vectorField& Foam::ggiPolyPatch::reconFaceCellCentres() const
{
    if (!reconFaceCellCentresPtr_)
    {
        calcReconFaceCellCentres();
    }

    return *reconFaceCellCentresPtr_;
}


void Foam::ggiPolyPatch::initAddressing()
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
            sendAddr();
        }
    }

    polyPatch::initAddressing();
}


void Foam::ggiPolyPatch::calcAddressing()
{
    polyPatch::calcAddressing();
}


void Foam::ggiPolyPatch::initGeometry()
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


void Foam::ggiPolyPatch::calcGeometry()
{
    polyPatch::calcGeometry();

    // Note: Calculation of transforms must be forced before the
    // reconFaceCellCentres in order to correctly set the transformation
    // in the interpolation routines.
    // HJ, 3/Jul/2009
}


void Foam::ggiPolyPatch::initMovePoints(const pointField& p)
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
            -separation()
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
            sendAddr();
        }
    }

    if (active() && master())
    {
        reconFaceCellCentres();
    }

    polyPatch::initMovePoints(p);
}


void Foam::ggiPolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);
}


void Foam::ggiPolyPatch::initUpdateMesh()
{
    polyPatch::initUpdateMesh();
}


void Foam::ggiPolyPatch::updateMesh()
{
    polyPatch::updateMesh();
    clearOut();
}


void Foam::ggiPolyPatch::calcTransforms() const
{
    // Simplest mixing plane: no transform or separation.  HJ, 24/Oct/2008
    forwardT_.setSize(0);
    reverseT_.setSize(0);
    separation_.setSize(0);

    if (debug > 1 && master())
    {
        if
        (
            !empty()
         && patchToPatch().uncoveredMasterFaces().size() > 0
        )
        {
            // Write uncovered master faces
            Info<< "Writing uncovered master faces for patch "
                << name() << " as VTK." << endl;

            const polyMesh& mesh = boundaryMesh().mesh();

            fileName fvPath(mesh.time().path()/"VTK");
            mkDir(fvPath);

            indirectPrimitivePatch::writeVTK
            (
                fvPath/fileName("uncoveredGgiFaces" + name()),
                IndirectList<face>
                (
                    localFaces(),
                    patchToPatch().uncoveredMasterFaces()
                ),
                localPoints()
            );
        }

        if
        (
            !shadow().empty()
         && patchToPatch().uncoveredSlaveFaces().size() > 0
        )
        {
            // Write uncovered master faces
            Info<< "Writing uncovered shadow faces for patch "
                << shadowName() << " as VTK." << endl;

            const polyMesh& mesh = boundaryMesh().mesh();

            fileName fvPath(mesh.time().path()/"VTK");
            mkDir(fvPath);
            Pout<< "shadow().localFaces(): " << shadow().localFaces().size()
                << " patchToPatch().uncoveredSlaveFaces().size(): "
                << patchToPatch().uncoveredSlaveFaces().size()
                << " shadow().localPoints(): " << shadow().localPoints().size()
                << endl;
            indirectPrimitivePatch::writeVTK
            (
                fvPath/fileName("uncoveredGgiFaces" + shadowName()),
                IndirectList<face>
                (
                    shadow().localFaces(),
                    patchToPatch().uncoveredSlaveFaces()
                ),
                shadow().localPoints()
            );
        }

        // Check for bridge overlap
        if (!bridgeOverlap())
        {
            if
            (
                (
                    patchToPatch().uncoveredMasterFaces().size() > 0
                    && !empty()
                )
             || (
                    !shadow().empty()
                 && patchToPatch().uncoveredSlaveFaces().size() > 0
                )
            )
            {
                FatalErrorIn("label ggiPolyPatch::calcTransforms() const")
                    << "ggi patch " << name() << " with shadow "
                    << shadowName() << " has "
                    << patchToPatch().uncoveredMasterFaces().size()
                    << " uncovered master faces and "
                    << patchToPatch().uncoveredSlaveFaces().size()
                    << " uncovered slave faces.  Bridging is switched off. "
                    << abort(FatalError);
            }
        }
        else
        {
            InfoIn("label ggiPolyPatch::calcTransforms() const")
                << "ggi patch " << name() << " with shadow "
                << shadowName() << " has "
                << patchToPatch().uncoveredMasterFaces().size()
                << " uncovered master faces and "
                << patchToPatch().uncoveredSlaveFaces().size()
                << " uncovered slave faces.  Bridging is switched on. "
                << endl;
        }
    }
}


void Foam::ggiPolyPatch::initOrder(const primitivePatch&) const
{}


bool Foam::ggiPolyPatch::order
(
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size(), -1);
    rotation.setSize(pp.size(), 0);

    // Nothing changes
    return false;
}


void Foam::ggiPolyPatch::syncOrder() const
{}


void Foam::ggiPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("shadowPatch") << shadowName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("zone") << zoneName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("bridgeOverlap") << bridgeOverlap_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
