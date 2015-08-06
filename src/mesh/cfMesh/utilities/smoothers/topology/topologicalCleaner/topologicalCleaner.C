/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

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
