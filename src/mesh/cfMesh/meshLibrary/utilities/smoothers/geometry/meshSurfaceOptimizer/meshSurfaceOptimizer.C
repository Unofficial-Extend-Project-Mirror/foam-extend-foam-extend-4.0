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
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "meshSurfaceMapper.H"
#include "polyMeshGenChecks.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceOptimizer::classifySurfaceVertices()
{
    const labelHashSet& corners = partitionerPtr_->corners();
    const labelHashSet& edgePoints = partitionerPtr_->edgePoints();

    //- set all vertices to partition
    vertexType_ = PARTITION;

    //- set corners
    forAllConstIter(labelHashSet, corners, it)
        vertexType_[it.key()] = CORNER;

    //- set edges
    forAllConstIter(labelHashSet, edgePoints, it)
        vertexType_[it.key()] = EDGE;

    if( Pstream::parRun() )
    {
        //- mark nodes at parallel boundaries
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();

        forAllConstIter(Map<label>, globalToLocal, iter)
        {
            const label bpI = iter();

            vertexType_[bpI] |= PROCBND;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceOptimizer::meshSurfaceOptimizer(const meshSurfaceEngine& surface)
:
    surfaceEngine_(surface),
    vertexType_(surface.boundaryPoints().size()),
    partitionerPtr_(new meshSurfacePartitioner(surface)),
    deletePartitioner_(true),
    octreePtr_(NULL),
    triMeshPtr_(NULL),
    enforceConstraints_(false),
    badPointsSubsetName_("invertedBoundaryPoints")
{
    classifySurfaceVertices();
}

meshSurfaceOptimizer::meshSurfaceOptimizer(const meshSurfacePartitioner& mPart)
:
    surfaceEngine_(mPart.surfaceEngine()),
    vertexType_(surfaceEngine_.boundaryPoints().size()),
    partitionerPtr_(&mPart),
    deletePartitioner_(true),
    octreePtr_(NULL),
    triMeshPtr_(NULL),
    enforceConstraints_(false),
    badPointsSubsetName_("invertedBoundaryPoints")
{
    classifySurfaceVertices();
}

meshSurfaceOptimizer::meshSurfaceOptimizer
(
    const meshSurfaceEngine& surface,
    const meshOctree& octree
)
:
    surfaceEngine_(surface),
    vertexType_(surface.boundaryPoints().size()),
    partitionerPtr_(new meshSurfacePartitioner(surface)),
    deletePartitioner_(true),
    octreePtr_(&octree),
    triMeshPtr_(NULL),
    enforceConstraints_(false),
    badPointsSubsetName_("invertedBoundaryPoints")
{
    classifySurfaceVertices();
}

meshSurfaceOptimizer::meshSurfaceOptimizer
(
    const meshSurfacePartitioner& partitioner,
    const meshOctree& octree
)
:
    surfaceEngine_(partitioner.surfaceEngine()),
    vertexType_(surfaceEngine_.boundaryPoints().size()),
    partitionerPtr_(&partitioner),
    deletePartitioner_(false),
    octreePtr_(&octree),
    triMeshPtr_(NULL),
    enforceConstraints_(false),
    badPointsSubsetName_("invertedBoundaryPoints")
{
    classifySurfaceVertices();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceOptimizer::~meshSurfaceOptimizer()
{
    deleteDemandDrivenData(triMeshPtr_);

    if( deletePartitioner_ )
        deleteDemandDrivenData(partitionerPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceOptimizer::removeUserConstraints()
{
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(vertexType_, bpI)
        if( vertexType_[bpI] & LOCKED )
            vertexType_[bpI] ^= LOCKED;
}

void meshSurfaceOptimizer::enforceConstraints(const word subsetName)
{
    enforceConstraints_ = true;

    badPointsSubsetName_ = subsetName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
