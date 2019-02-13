/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "demandDrivenData.H"
#include "topologicalCleaner.H"
#include "decomposeCells.H"
#include "checkCellConnectionsOverFaces.H"
#include "meshSurfaceEngine.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void topologicalCleaner::decomposeCells()
{
    if( !changed_ )
         return;

    Foam::decomposeCells dc(mesh_);
    dc.decomposeMesh(decomposeCell_);

}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from points, cells, boundary faces, and octree
topologicalCleaner::topologicalCleaner
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
    changed_(false),
    decomposeCell_(mesh.cells().size(), false)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

topologicalCleaner::~topologicalCleaner()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool topologicalCleaner::cleanTopology()
{
    checkInvalidConnectionsForVerticesCells();

    checkInvalidConnectionsForVerticesFaces();

    checkNonConsecutiveBoundaryVertices();

    checkNonMappableCells();

    checkNonMappableFaces();

    decomposeCells();

    if( checkCellConnectionsOverFaces(mesh_).checkCellGroups() )
        changed_ = true;

    return changed_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
