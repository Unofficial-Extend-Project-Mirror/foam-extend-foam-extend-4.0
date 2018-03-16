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

#include "fvMesh.H"
#include "fvBoundaryMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void fvBoundaryMesh::addFvPatches()
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
    
    setSize(bMesh.size());

    // Set boundary patches, using the patches added to the polyMesh
    // Bug fix.  HJ, 1/Mar/2018
    fvPatchList& Patches = *this;

    forAll(Patches, patchI)
    {
        Patches.set(patchI, fvPatch::New(bMesh[patchI], *this));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fvBoundaryMesh::fvBoundaryMesh
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

void fvBoundaryMesh::movePoints()
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


lduInterfacePtrsList fvBoundaryMesh::interfaces() const
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


void fvBoundaryMesh::readUpdate()
{
    clear();
    addFvPatches();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
