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

#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceEngine::meshSurfaceEngine(polyMeshGen& mesh)
:
    mesh_(mesh),
    activePatch_(-1),
    boundaryPointsPtr_(nullptr),
    boundaryFacesPtr_(nullptr),
    boundaryFacePatchPtr_(nullptr),
    boundaryFaceOwnersPtr_(nullptr),
    pointFacesPtr_(nullptr),
    pointInFacePtr_(nullptr),
    pointPatchesPtr_(nullptr),
    bppPtr_(nullptr),
    pointPointsPtr_(nullptr),
    edgesPtr_(nullptr),
    bpEdgesPtr_(nullptr),
    edgeFacesPtr_(nullptr),
    faceEdgesPtr_(nullptr),
    edgePatchesPtr_(nullptr),
    faceFacesPtr_(nullptr),
    pointNormalsPtr_(nullptr),
    faceNormalsPtr_(nullptr),
    faceCentresPtr_(nullptr),

    globalBoundaryPointLabelPtr_(nullptr),
    globalBoundaryPointToLocalPtr_(nullptr),
    bpProcsPtr_(nullptr),
    bpNeiProcsPtr_(nullptr),
    globalBoundaryEdgeLabelPtr_(nullptr),
    globalBoundaryEdgeToLocalPtr_(nullptr),
    beProcsPtr_(nullptr),
    beNeiProcsPtr_(nullptr),
    otherEdgeFaceAtProcPtr_(nullptr),
    otherEdgeFacePatchPtr_(nullptr),
    globalBoundaryFaceLabelPtr_(nullptr)
{
    calculateBoundaryFaces();
    calculateBoundaryNodes();
}

meshSurfaceEngine::meshSurfaceEngine(polyMeshGen &mesh, const label patchI)
:
    mesh_(mesh),
    activePatch_(patchI),
    boundaryPointsPtr_(nullptr),
    boundaryFacesPtr_(nullptr),
    boundaryFacePatchPtr_(nullptr),
    boundaryFaceOwnersPtr_(nullptr),
    pointFacesPtr_(nullptr),
    pointInFacePtr_(nullptr),
    pointPatchesPtr_(nullptr),
    bppPtr_(nullptr),
    pointPointsPtr_(nullptr),
    edgesPtr_(nullptr),
    bpEdgesPtr_(nullptr),
    edgeFacesPtr_(nullptr),
    faceEdgesPtr_(nullptr),
    edgePatchesPtr_(nullptr),
    faceFacesPtr_(nullptr),
    pointNormalsPtr_(nullptr),
    faceNormalsPtr_(nullptr),
    faceCentresPtr_(nullptr),

    globalBoundaryPointLabelPtr_(nullptr),
    globalBoundaryPointToLocalPtr_(nullptr),
    bpProcsPtr_(nullptr),
    bpNeiProcsPtr_(nullptr),
    globalBoundaryEdgeLabelPtr_(nullptr),
    globalBoundaryEdgeToLocalPtr_(nullptr),
    beProcsPtr_(nullptr),
    beNeiProcsPtr_(nullptr),
    otherEdgeFaceAtProcPtr_(nullptr),
    otherEdgeFacePatchPtr_(nullptr),
    globalBoundaryFaceLabelPtr_(nullptr)
{
    calculateBoundaryFaces();
    calculateBoundaryNodes();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceEngine::~meshSurfaceEngine()
{
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceEngine::clearOut()
{
    deleteDemandDrivenData(boundaryPointsPtr_);
    deleteDemandDrivenData(boundaryFacesPtr_);
    deleteDemandDrivenData(boundaryFacePatchPtr_);
    deleteDemandDrivenData(boundaryFaceOwnersPtr_);
    deleteDemandDrivenData(pointFacesPtr_);
    deleteDemandDrivenData(pointInFacePtr_);
    deleteDemandDrivenData(pointPatchesPtr_);
    deleteDemandDrivenData(bppPtr_);
    deleteDemandDrivenData(pointPointsPtr_);
    deleteDemandDrivenData(pointNormalsPtr_);
    deleteDemandDrivenData(faceNormalsPtr_);
    deleteDemandDrivenData(faceCentresPtr_);
    deleteDemandDrivenData(edgesPtr_);
    deleteDemandDrivenData(bpEdgesPtr_);
    deleteDemandDrivenData(edgeFacesPtr_);
    deleteDemandDrivenData(faceEdgesPtr_);
    deleteDemandDrivenData(edgePatchesPtr_);
    deleteDemandDrivenData(faceFacesPtr_);

    deleteDemandDrivenData(globalBoundaryPointLabelPtr_);
    deleteDemandDrivenData(globalBoundaryPointToLocalPtr_);
    deleteDemandDrivenData(bpProcsPtr_);
    deleteDemandDrivenData(bpNeiProcsPtr_);
    deleteDemandDrivenData(globalBoundaryEdgeLabelPtr_);
    deleteDemandDrivenData(globalBoundaryEdgeToLocalPtr_);
    deleteDemandDrivenData(beProcsPtr_);
    deleteDemandDrivenData(beNeiProcsPtr_);
    deleteDemandDrivenData(otherEdgeFaceAtProcPtr_);
    deleteDemandDrivenData(otherEdgeFacePatchPtr_);
    deleteDemandDrivenData(globalBoundaryFaceLabelPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
