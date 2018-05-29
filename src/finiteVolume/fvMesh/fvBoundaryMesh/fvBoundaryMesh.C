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

#include "fvMesh.H"
#include "fvBoundaryMesh.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvBoundaryMesh::addFvPatches()
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();

    // Clear existing patches and resize the list
    clear();
    setSize(bMesh.size());

    // Set boundary patches, using the patches added to the polyMesh
    // Bug fix.  HJ, 1/Mar/2018
    fvPatchList& Patches = *this;

    forAll (Patches, patchI)
    {
        Patches.set(patchI, fvPatch::New(bMesh[patchI], *this));
    }
}


void Foam::fvBoundaryMesh::resetFvPatches(const boolList& resetFvPatchFlag)
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();

    if (resetFvPatchFlag.size() != bMesh.size())
    {
        FatalErrorIn("void resetFvPatches(const boolList& resetFvPatchFlag)")
            << "Incorrect size of reset list.  Boundary size: "
            << bMesh.size() << " reset size: " << resetFvPatchFlag.size()
            << abort(FatalError);
    }

    // Reset list size.  This will delete pointers to patches
    // if the list is truncated
    setSize(bMesh.size());

    // Set boundary patches, using the patches added to the polyMesh
    // Bug fix.  HJ, 1/Mar/2018
    fvPatchList& Patches = *this;

    // In order to preserve patch links on resize of boundary,
    // only reset the empty slots
    forAll (Patches, patchI)
    {
        if (resetFvPatchFlag[patchI])
        {
            // Set new patch.  This also deletes old pointer
            Patches.set(patchI, fvPatch::New(bMesh[patchI], *this));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvBoundaryMesh::fvBoundaryMesh
(
    const fvMesh& m
)
:
    fvPatchList(m.boundaryMesh().size()),
    mesh_(m)
{
    addFvPatches();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvBoundaryMesh::movePoints()
{
    forAll(*this, patchi)
    {
        operator[](patchi).initMovePoints();
    }

    forAll(*this, patchi)
    {
        operator[](patchi).movePoints();
    }
}


Foam::lduInterfacePtrsList Foam::fvBoundaryMesh::interfaces() const
{
    lduInterfacePtrsList interfaces(size());

    forAll (interfaces, patchi)
    {
        if (isA<lduInterface>(this->operator[](patchi)))
        {
            interfaces.set
            (
                patchi,
               &refCast<const lduInterface>(this->operator[](patchi))
            );
        }
    }

    return interfaces;
}


void Foam::fvBoundaryMesh::readUpdate()
{
    clear();
    addFvPatches();
}


// ************************************************************************* //
