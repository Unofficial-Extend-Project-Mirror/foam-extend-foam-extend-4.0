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

label meshSurfaceOptimizer::findBadFaces
(
    labelHashSet& badFaces,
    boolList& changedFace
) const
{
    badFaces.clear();

    const polyMeshGen& mesh = surfaceEngine_.mesh();

    polyMeshGenChecks::checkFacePyramids
    (
        mesh,
        false,
        VSMALL,
        &badFaces,
        &changedFace
    );

    polyMeshGenChecks::checkCellPartTetrahedra
    (
        mesh,
        false,
        VSMALL,
        &badFaces,
        &changedFace
    );

    polyMeshGenChecks::checkFaceAreas
    (
        mesh,
        false,
        VSMALL,
        &badFaces,
        &changedFace
    );

    const label nBadFaces = returnReduce(badFaces.size(), sumOp<label>());

    return nBadFaces;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceOptimizer::meshSurfaceOptimizer
(
    const meshSurfaceEngine& surface,
    const meshOctree& octree
)
:
    surfaceEngine_(surface),
    meshOctree_(octree),
    vertexType_(surface.boundaryPoints().size()),
    partitionerPtr_(new meshSurfacePartitioner(surface)),
    deletePartitioner_(true),
    triMeshPtr_(NULL),
    enforceConstraints_(false),
    badPointsSubsetName_()
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
    meshOctree_(octree),
    vertexType_(surfaceEngine_.boundaryPoints().size()),
    partitionerPtr_(&partitioner),
    deletePartitioner_(false),
    triMeshPtr_(NULL),
    enforceConstraints_(false),
    badPointsSubsetName_()
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

void meshSurfaceOptimizer::enforceConstraints(const word subsetName)
{
    enforceConstraints_ = true;

    badPointsSubsetName_ = subsetName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
