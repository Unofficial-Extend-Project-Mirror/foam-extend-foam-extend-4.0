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

#include "meshSurfaceMapper.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "triSurf.H"
#include "triSurfacePartitioner.H"
#include "demandDrivenData.H"
#include "meshOctree.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper::createMeshSurfacePartitioner() const
{
    surfaceEnginePartitionerPtr_ = new meshSurfacePartitioner(surfaceEngine_);
}

void meshSurfaceMapper::createTriSurfacePartitioner() const
{
    surfPartitionerPtr_ = new triSurfacePartitioner(meshOctree_.surface());
}
    
void meshSurfaceMapper::clearOut()
{
    if( deletePartitioner_ )
        deleteDemandDrivenData(surfaceEnginePartitionerPtr_);
    deleteDemandDrivenData(surfPartitionerPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceMapper::meshSurfaceMapper
(
    const meshSurfaceEngine& mse,
    const meshOctree& octree
)
:
    surfaceEngine_(mse),
    meshOctree_(octree),
    surfaceEnginePartitionerPtr_(NULL),
    deletePartitioner_(true),
    surfPartitionerPtr_(NULL)
{
    if( Pstream::parRun() )
    {
        //- allocate bpAtProcs and other addressing
        //- this is done here to prevent possible deadlocks
        surfaceEngine_.bpAtProcs();
    }
}

meshSurfaceMapper::meshSurfaceMapper
(
    const meshSurfacePartitioner& mPart,
    const meshOctree& octree
)
:
    surfaceEngine_(mPart.surfaceEngine()),
    meshOctree_(octree),
    surfaceEnginePartitionerPtr_(&mPart),
    deletePartitioner_(false),
    surfPartitionerPtr_(NULL)
{
    if( Pstream::parRun() )
    {
        //- allocate bpAtProcs and other addressing
        //- this is done here to prevent possible deadlocks
        surfaceEngine_.bpAtProcs();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceMapper::~meshSurfaceMapper()
{
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
