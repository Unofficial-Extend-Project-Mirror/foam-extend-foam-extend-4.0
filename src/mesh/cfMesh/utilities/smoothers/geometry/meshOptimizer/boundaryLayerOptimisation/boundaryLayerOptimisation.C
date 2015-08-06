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

#include "boundaryLayerOptimisation.H"
#include "meshSurfacePartitioner.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

boundaryLayerOptimisation::boundaryLayerOptimisation(polyMeshGen& mesh)
:
    mesh_(mesh),
    meshSurfacePtr_(new meshSurfaceEngine(mesh)),
    deleteMeshSurface_(true),
    partitionerPtr_(NULL),
    hairEdges_(),
    hairEdgesAtBndPoint_(),
    hairEdgesNearHairEdge_(),
    isBndLayerBase_(),
    isExitFace_(),
    hairEdgeType_(),
    thinnedHairEdge_(),
    maxNumIterations_(5),
    nSmoothNormals_(5),
    relThicknessTol_(0.1),
    featureSizeFactor_(0.3),
    reCalculateNormals_(true)
{
    calculateHairEdges();
}

boundaryLayerOptimisation::boundaryLayerOptimisation
(
    polyMeshGen& mesh,
    const meshSurfaceEngine& mse
)
:
    mesh_(mesh),
    meshSurfacePtr_(&mse),
    deleteMeshSurface_(false),
    partitionerPtr_(NULL),
    hairEdges_(),
    hairEdgesAtBndPoint_(),
    hairEdgesNearHairEdge_(),
    isBndLayerBase_(),
    isExitFace_(),
    hairEdgeType_(),
    thinnedHairEdge_(),
    maxNumIterations_(5),
    nSmoothNormals_(5),
    relThicknessTol_(0.15),
    featureSizeFactor_(0.3),
    reCalculateNormals_(true)
{
    calculateHairEdges();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

boundaryLayerOptimisation::~boundaryLayerOptimisation()
{
    deleteDemandDrivenData(partitionerPtr_);

    if( deleteMeshSurface_ )
        deleteDemandDrivenData(meshSurfacePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayerOptimisation::setMaxNumIterations
(
    const label maxNumIterations
)
{
    maxNumIterations_ = maxNumIterations;
}

void boundaryLayerOptimisation::setNumNormalsSmoothingIterations
(
    const label nSmoothNormals
)
{
    nSmoothNormals_ = nSmoothNormals;
}

void boundaryLayerOptimisation::recalculateNormals(const bool shallRecalculate)
{
    reCalculateNormals_ = shallRecalculate;
}

void boundaryLayerOptimisation::setRelativeThicknessTolerance
(
    const scalar relThicknessTol
)
{
    relThicknessTol_ = relThicknessTol;
}

void boundaryLayerOptimisation::setFeatureSizeFactor
(
    const scalar featureSizeFactor
)
{
    featureSizeFactor_ = featureSizeFactor;
}

const edgeLongList& boundaryLayerOptimisation::hairEdges() const
{
    return hairEdges_;
}

const VRWGraph& boundaryLayerOptimisation::hairEdgesAtBndPoint() const
{
    return hairEdgesAtBndPoint_;
}

const boolList& boundaryLayerOptimisation::isBaseFace() const
{
    return isBndLayerBase_;
}

const boolList& boundaryLayerOptimisation::isExitFace() const
{
    return isExitFace_;
}

void boundaryLayerOptimisation::readSettings
(
    const dictionary& meshDict,
    boundaryLayerOptimisation& blOptimisation
)
{
    if( meshDict.found("boundaryLayers") )
    {
        const dictionary& layersDict = meshDict.subDict("boundaryLayers");

        if( layersDict.found("optimiseLayer") )
        {
            const bool smoothLayers =
                readBool(layersDict.lookup("optimiseLayer"));

            if( !smoothLayers )
                return;
        }

        if( layersDict.found("optimisationParameters") )
        {
            const dictionary& optParams =
                layersDict.subDict("optimisationParameters");

            if( optParams.found("recalculateNormals") )
            {
                const bool recalculateNormals =
                    readBool(optParams.lookup("recalculateNormals"));

                blOptimisation.recalculateNormals(recalculateNormals);
            }

            if( optParams.found("nSmoothNormals") )
            {
                const label nSmoothNormals =
                    readLabel(optParams.lookup("nSmoothNormals"));

                blOptimisation.setNumNormalsSmoothingIterations(nSmoothNormals);
            }

            if( optParams.found("featureSizeFactor") )
            {
                const scalar featureSizeFactor =
                    readScalar(optParams.lookup("featureSizeFactor"));

                if( featureSizeFactor >= 1.0 || featureSizeFactor < 0.0 )
                    FatalErrorIn
                    (
                        "void boundaryLayerOptimisation::optimiseLayer"
                        "(const dictionary&, boundaryLayerOptimisation&)"
                    ) << "Feature size factor is out"
                      << " of a valid range 0 to 1" << exit(FatalError);

                blOptimisation.setFeatureSizeFactor(featureSizeFactor);
            }

            if( optParams.found("relThicknessTol") )
            {
                const scalar relThicknessTol =
                    readScalar(optParams.lookup("relThicknessTol"));

                if( relThicknessTol >= 1.0 || relThicknessTol < 0.0 )
                    FatalErrorIn
                    (
                        "void boundaryLayerOptimisation::optimiseLayer"
                        "(const dictionary&, boundaryLayerOptimisation&)"
                    ) << "Relative thickness tolerance is out"
                      << " of a valid range 0 to 1" << exit(FatalError);

                blOptimisation.setRelativeThicknessTolerance(relThicknessTol);
            }

            if( optParams.found("maxNumIterations") )
            {
                const label maxNumIterations =
                    readLabel(optParams.lookup("maxNumIterations"));

                blOptimisation.setMaxNumIterations(maxNumIterations);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
