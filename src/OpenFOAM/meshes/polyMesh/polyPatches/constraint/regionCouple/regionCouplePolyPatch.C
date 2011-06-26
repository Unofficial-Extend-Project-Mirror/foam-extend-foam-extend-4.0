/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-7 H. Jasak All rights reserved
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
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "regionCouplePolyPatch.H"
#include "polyBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "demandDrivenData.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyPatchID.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCouplePolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, regionCouplePolyPatch, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::regionCouplePolyPatch::calcInterpolation() const
{
    // Create patch-to-patch interpolation
    if (patchToPatchPtr_)
    {
        FatalErrorIn("void regionCouplePolyPatch::calcInterpolation() const")
            << "Patch to patch interpolation already calculated"
            << abort(FatalError);
    }

    // Get shadow region
    const polyMesh& sr = shadowRegion();

    // Grab shadow patch index
    polyPatchID shadow(shadowPatchName_, sr.boundaryMesh());

    if (!shadow.active())
    {
        FatalErrorIn("void regionCouplePolyPatch::calcInterpolation() const")
            << "Shadow patch name " << shadowPatchName_
            << " not found.  Please check your regionCouple interface."
            << abort(FatalError);
    }

    shadowIndex_ = shadow.index();

    // Check the other side is a regionCouple
    if (!isType<regionCouplePolyPatch>(sr.boundaryMesh()[shadowIndex_]))
    {
        FatalErrorIn("void regionCouplePolyPatch::calcInterpolation() const")
            << "Shadow of regionCouple patch " << name()
            << " named " << shadowPatchName()
            << " is not a regionCouple." << nl
            << "This is not allowed.  Please check your mesh definition."
            << abort(FatalError);
    }

    patchToPatchPtr_ =
        new patchToPatchInterpolation
        (
            sr.boundaryMesh()[shadowIndex_],
            *this,
            intersection::VISIBLE
        );
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


const Foam::patchToPatchInterpolation&
Foam::regionCouplePolyPatch::patchToPatch() const
{
    if (!attached_)
    {
        FatalErrorIn
        (
            "const patchToPatchInterpolation& "
            "regionCouplePolyPatch::patchToPatch() const"
        )   << "Requesting patchToPatchInterpolation in detached state"
            << abort(FatalError);
    }

    if (!patchToPatchPtr_)
    {
        calcInterpolation();
    }

    return *patchToPatchPtr_;
}


void Foam::regionCouplePolyPatch::calcReconFaceCellCentres() const
{
    // Create neighbouring face centres using interpolation

    // Get shadow region
    const polyMesh& sr = shadowRegion();

    const label shadowID = shadowIndex();

    // Reconstruct the shadow cell face centres
    vectorField localCtrs = faceCellCentres();
    vectorField reconCtrs =
        patchToPatch().faceInterpolate
        (
            sr.boundaryMesh()[shadowID].faceCellCentres()
        );

    // Calculate reconstructed centres by eliminating non-orthogonality
    const vectorField& n = faceNormals();

    reconFaceCellCentresPtr_ =
        new vectorField(localCtrs + n*(n & (reconCtrs - localCtrs)));
}


void Foam::regionCouplePolyPatch::clearGeom() const
{
    deleteDemandDrivenData(reconFaceCellCentresPtr_);
}


void Foam::regionCouplePolyPatch::clearOut() const
{
    clearGeom();

    deleteDemandDrivenData(patchToPatchPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::regionCouplePolyPatch::regionCouplePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& shadowRegionName,
    const word& shadowPatchName,
    const bool attached,
    const bool isWall
)
:
    coupledPolyPatch(name, size, start, index, bm),
    shadowRegionName_(shadowRegionName),
    shadowPatchName_(shadowPatchName),
    attached_(attached),
    isWall_(isWall),
    shadowIndex_(-1),
    patchToPatchPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
{}


// Construct from dictionary
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
    attached_(dict.lookup("attached")),
    isWall_(dict.lookup("isWall")),
    shadowIndex_(-1),
    patchToPatchPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
{}


//- Construct as copy, resetting the boundary mesh
Foam::regionCouplePolyPatch::regionCouplePolyPatch
(
    const regionCouplePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    shadowRegionName_(pp.shadowRegionName_),
    shadowPatchName_(pp.shadowPatchName_),
    attached_(pp.attached_),
    isWall_(pp.isWall_),
    shadowIndex_(-1),
    patchToPatchPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
{}


//- Construct as copy, resetting the face list and boundary mesh data
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
    attached_(pp.attached_),
    isWall_(pp.isWall_),
    shadowIndex_(-1),
    patchToPatchPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
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
    return boundaryMesh().mesh().db().parent().
        objectRegistry::lookupObject<polyMesh>
        (
            shadowRegionName()
        );
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
        clearGeom();
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
        clearGeom();
    }
}


Foam::label Foam::regionCouplePolyPatch::shadowIndex() const
{
    if (shadowIndex_ == -1)
    {
        calcInterpolation();
    }

    return shadowIndex_;
}


const Foam::regionCouplePolyPatch& Foam::regionCouplePolyPatch::shadow() const
{
    return refCast<const regionCouplePolyPatch>
    (
        shadowRegion().boundaryMesh()[shadowIndex()]
    );
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
    polyPatch::initAddressing();
}


void Foam::regionCouplePolyPatch::calcAddressing()
{
    polyPatch::calcAddressing();
}


void Foam::regionCouplePolyPatch::initGeometry()
{
    polyPatch::initGeometry();
}


void Foam::regionCouplePolyPatch::calcGeometry()
{
    polyPatch::calcGeometry();
    // Reconstruct the cell face centres
//     reconFaceCellCentres();
}


void Foam::regionCouplePolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::initMovePoints(p);
}


void Foam::regionCouplePolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);

    // Clear reconstructed face centres
    deleteDemandDrivenData(reconFaceCellCentresPtr_);

    if (patchToPatchPtr_)
    {
        patchToPatchPtr_->movePoints();
    }
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


void Foam::regionCouplePolyPatch::calcTransformTensors
(
    const vectorField& Cf,
    const vectorField& Cr,
    const vectorField& nf,
    const vectorField& nry
) const
{
    FatalErrorIn("void regionCouplePolyPatch::calcTransformTensors")
        << "Not ready"
        << abort(FatalError);
}


void Foam::regionCouplePolyPatch::initOrder(const primitivePatch& pp) const
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
    os.writeKeyword("attached")
        << attached_ << token::END_STATEMENT << nl;
    os.writeKeyword("isWall")
        << isWall_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
