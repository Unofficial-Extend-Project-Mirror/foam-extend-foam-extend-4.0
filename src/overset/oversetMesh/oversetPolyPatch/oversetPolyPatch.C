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

#include "oversetPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "foamTime.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, oversetPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        oversetPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::oversetPolyPatch::initAddressing()
{
    polyPatch::initAddressing();
}


void Foam::oversetPolyPatch::calcAddressing()
{
    polyPatch::calcAddressing();
}


void Foam::oversetPolyPatch::initGeometry()
{
    polyPatch::initGeometry();
}


void Foam::oversetPolyPatch::calcGeometry()
{
    polyPatch::calcGeometry();
}


void Foam::oversetPolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::movePoints(p);
    oversetPolyPatch::initGeometry();
}


void Foam::oversetPolyPatch::movePoints(const pointField&)
{
    oversetPolyPatch::calcGeometry();
}


void Foam::oversetPolyPatch::initUpdateMesh()
{
    polyPatch::initUpdateMesh();
}


void Foam::oversetPolyPatch::updateMesh()
{
    polyPatch::updateMesh();
}


void Foam::oversetPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetPolyPatch::oversetPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, size, start, index, bm)
{}


Foam::oversetPolyPatch::oversetPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm)
{}


Foam::oversetPolyPatch::oversetPolyPatch
(
    const oversetPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm)
{}


Foam::oversetPolyPatch::oversetPolyPatch
(
    const oversetPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::oversetPolyPatch::~oversetPolyPatch()
{
//     clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::oversetPolyPatch::initOrder(const primitivePatch& pp) const
{}


// Return new ordering. Ordering is -faceMap: for every face index
// the new face -rotation:for every new face the clockwise shift
// of the original face. Return false if nothing changes (faceMap
// is identity, rotation is 0)
bool Foam::oversetPolyPatch::order
(
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    return false;
}


void Foam::oversetPolyPatch::syncOrder() const
{}


// ************************************************************************* //
