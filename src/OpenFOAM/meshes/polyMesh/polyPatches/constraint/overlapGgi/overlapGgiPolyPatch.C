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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.

\*---------------------------------------------------------------------------*/

#include "overlapGgiPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "demandDrivenData.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"
#include "SubField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(overlapGgiPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, overlapGgiPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, overlapGgiPolyPatch, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::overlapGgiPolyPatch::active() const
{
    polyPatchID shadow(shadowName_, boundaryMesh());
    faceZoneID zone(zoneName_, boundaryMesh().mesh().faceZones());

    return shadow.active() && zone.active();
}


void Foam::overlapGgiPolyPatch::calcLocalParallel() const
{
    // Calculate patch-to-zone addressing
    if (localParallelPtr_)
    {
        FatalErrorIn("void overlapGgiPolyPatch::calcLocalParallel() const")
            << "Local parallel switch already calculated"
            << abort(FatalError);
    }

    // If running in serial, all GGIs are local parallel
    // HJ, 1/Jun/2011
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
        Info<< "Overlap GGI patch Master: " << name()
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::overlapGgiPolyPatch::overlapGgiPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, size, start, index, bm),
    shadowName_(word::null),
    zoneName_(word::null),
    shadowIndex_(-1),
    zoneIndex_(-1),
    rotationAxis_(vector(0, 0, 1)),
    nCopies_(0),
    expandedMasterPtr_(NULL),
    expandedSlavePtr_(NULL),
    patchToPatchPtr_(NULL),
    localParallelPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
{}


// Construct from components
Foam::overlapGgiPolyPatch::overlapGgiPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& shadowName,
    const word& zoneName,
    const vector& axis,
    const scalar nCopies
)
:
    coupledPolyPatch(name, size, start, index, bm),
    shadowName_(shadowName),
    zoneName_(zoneName),
    shadowIndex_(-1),
    zoneIndex_(-1),
    rotationAxis_(axis),
    nCopies_(nCopies),
    expandedMasterPtr_(NULL),
    expandedSlavePtr_(NULL),
    patchToPatchPtr_(NULL),
    localParallelPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
{}


// Construct from dictionary
Foam::overlapGgiPolyPatch::overlapGgiPolyPatch
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
    shadowIndex_(-1),
    zoneIndex_(-1),
    rotationAxis_(dict.lookup("rotationAxis")),
    nCopies_(readScalar(dict.lookup("nCopies"))),
    expandedMasterPtr_(NULL),
    expandedSlavePtr_(NULL),
    patchToPatchPtr_(NULL),
    localParallelPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
{}


//- Construct as copy, resetting the boundary mesh
Foam::overlapGgiPolyPatch::overlapGgiPolyPatch
(
    const overlapGgiPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    shadowName_(pp.shadowName_),
    zoneName_(pp.zoneName_),
    shadowIndex_(-1),
    zoneIndex_(-1),
    rotationAxis_(pp.rotationAxis_),
    nCopies_(pp.nCopies_),
    expandedMasterPtr_(NULL),
    expandedSlavePtr_(NULL),
    patchToPatchPtr_(NULL),
    localParallelPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
{}


//- Construct as copy, resetting the face list and boundary mesh data
Foam::overlapGgiPolyPatch::overlapGgiPolyPatch
(
    const overlapGgiPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    shadowName_(pp.shadowName_),
    zoneName_(pp.zoneName_),
    shadowIndex_(-1),
    zoneIndex_(-1),
    rotationAxis_(pp.rotationAxis_),
    nCopies_(pp.nCopies_),
    expandedMasterPtr_(NULL),
    expandedSlavePtr_(NULL),
    patchToPatchPtr_(NULL),
    localParallelPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::overlapGgiPolyPatch::~overlapGgiPolyPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::overlapGgiPolyPatch::shadowIndex() const
{
    if (shadowIndex_ == -1 && shadowName_ != word::null)
    {
        // Grab shadow patch index
        polyPatchID shadow(shadowName_, boundaryMesh());

        if (!shadow.active())
        {
            FatalErrorIn("label overlapGgiPolyPatch::shadowIndex() const")
                << "Shadow patch name " << shadowName_
                << " not found.  Please check your GGI interface definition."
                << abort(FatalError);
        }

        shadowIndex_ = shadow.index();

        // Check the other side is a ggi
        if (!isA<overlapGgiPolyPatch>(boundaryMesh()[shadowIndex_]))
        {
            FatalErrorIn("label overlapGgiPolyPatch::shadowIndex() const")
                << "Shadow of ggi patch " << name()
                << " named " << shadowName() << " is not a ggi.  Type: "
                << boundaryMesh()[shadowIndex_].type() << nl
                << "This is not allowed.  Please check your mesh definition."
                << abort(FatalError);
        }

        if (index() == shadowIndex_)
        {
            FatalErrorIn("label overlapGgiPolyPatch::shadowIndex() const")
                << "ggi patch " << name() << " created as its own shadow"
                << abort(FatalError);
        }
    }

    return shadowIndex_;
}

Foam::label Foam::overlapGgiPolyPatch::zoneIndex() const
{
    if (zoneIndex_ == -1 && zoneName_ != Foam::word::null)
    {
        // Grab zone patch index
        faceZoneID zone(zoneName_, boundaryMesh().mesh().faceZones());

        if (!zone.active())
        {
            FatalErrorIn("label overlapGgiPolyPatch::zoneIndex() const")
                << "Face zone name " << zoneName_
                << " not found.  Please check your GGI interface definition."
                << abort(FatalError);
        }

        zoneIndex_ = zone.index();
    }

    return zoneIndex_;
}

const Foam::overlapGgiPolyPatch&
Foam::overlapGgiPolyPatch::shadow() const
{
    return refCast<const overlapGgiPolyPatch>(boundaryMesh()[shadowIndex()]);
}


const Foam::faceZone& Foam::overlapGgiPolyPatch::zone() const
{
    return boundaryMesh().mesh().faceZones()[zoneIndex()];
}


bool Foam::overlapGgiPolyPatch::localParallel() const
{
    // Calculate patch-to-zone addressing
    if (!localParallelPtr_)
    {
        calcLocalParallel();
    }

    return *localParallelPtr_;
}


const Foam::label& Foam::overlapGgiPolyPatch::nCopies() const
{
    // Read the number of copies to be made from the dictionary for the
    // expanded slave and expanded master to cover 360 degrees
    return nCopies_;
}


bool Foam::overlapGgiPolyPatch::master() const
{
    // The first overlapggi interface is master,second one is slave
    if (angle() == shadow().angle())
    {
        return index() < shadowIndex();
    }

    // Master is the one with the larger angle
    return angle() > shadow().angle();
}


// Write
void Foam::overlapGgiPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("rotationAxis") << rotationAxis_
        << token::END_STATEMENT << nl;
    os.writeKeyword("nCopies") << nCopies_
        << token::END_STATEMENT << nl;
    os.writeKeyword("shadowPatch") << shadowName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("zone") << zoneName_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
