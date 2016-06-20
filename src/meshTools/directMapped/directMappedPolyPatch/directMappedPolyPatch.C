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

#include "directMappedPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directMappedPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, directMappedPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, directMappedPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm),
    directMappedPatchBase(static_cast<const polyPatch&>(*this))
{}


Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const directMappedPatchBase::sampleMode mode,
    const word& samplePatch,
    const vectorField& offset,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm),
    directMappedPatchBase
    (
        static_cast<const polyPatch&>(*this),
        sampleRegion,
        mode,
        samplePatch,
        offset
    )
{}


Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const directMappedPatchBase::sampleMode mode,
    const word& samplePatch,
    const vector& offset,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm),
    directMappedPatchBase
    (
        static_cast<const polyPatch&>(*this),
        sampleRegion,
        mode,
        samplePatch,
        offset
    )
{}


Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, dict, index, bm),
    directMappedPatchBase(*this, dict)
{}


Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const directMappedPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    directMappedPatchBase(*this, pp)
{}


Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const directMappedPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    directMappedPatchBase(*this, pp)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::directMappedPolyPatch::~directMappedPolyPatch()
{
    directMappedPatchBase::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Initialise the calculation of the patch geometry
void Foam::directMappedPolyPatch::initGeometry()
{
    polyPatch::initGeometry();
}

//- Calculate the patch geometry
void Foam::directMappedPolyPatch::calcGeometry()
{
    polyPatch::calcGeometry();
}

//- Initialise the patches for moving points
void Foam::directMappedPolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::initMovePoints(p);

    // Force recalculation of mapping with new point position
    // Note: this uses parallel communications.  HJ, 13/Mar/2012
    directMappedPatchBase::clearOut();
    directMappedPatchBase::map();
}

//- Correct patches after moving points
void Foam::directMappedPolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);
}

//- Initialise the update of the patch topology
void Foam::directMappedPolyPatch::initUpdateMesh()
{
    polyPatch::initUpdateMesh();

    // Force recalculation of mapping with new point position
    // Note: this uses parallel communications.  HJ, 13/Mar/2012
    directMappedPatchBase::clearOut();

    // Only carry out mapping if the sampled region has been created already
    // DC, 04/Nov/2013
    if (boundaryMesh().mesh().time().found(sampleRegion()))
    {
        directMappedPatchBase::map();
    }
}

//- Update of the patch topology
void Foam::directMappedPolyPatch::updateMesh()
{
    polyPatch::updateMesh();
}


void Foam::directMappedPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    directMappedPatchBase::write(os);
}


// ************************************************************************* //
