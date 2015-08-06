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

#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const polyMeshGen& meshSurfaceEngine::mesh() const
{
    return mesh_;
}

inline const pointFieldPMG& meshSurfaceEngine::points() const
{
    return mesh_.points();
}

inline const faceListPMG& meshSurfaceEngine::faces() const
{
    return mesh_.faces();
}

inline const cellListPMG& meshSurfaceEngine::cells() const
{
    return mesh_.cells();
}

inline const labelList& meshSurfaceEngine::bp() const
{
    if( !bppPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelList& meshSurfaceEngine::bp() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateBoundaryFaces();
        calculateBoundaryNodes();
    }

    return *bppPtr_;
}

inline const labelList& meshSurfaceEngine::boundaryPoints() const
{
    if( !boundaryPointsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelList& meshSurfaceEngine::boundaryPoints() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateBoundaryNodes();
    }

    return *boundaryPointsPtr_;
}

inline const faceList::subList& meshSurfaceEngine::boundaryFaces() const
{
    if( !boundaryFacesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const faceList::subList&"
                "meshSurfaceEngine::boundaryFaces() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateBoundaryFaces();
    }

    return *boundaryFacesPtr_;
}

inline const labelList& meshSurfaceEngine::boundaryFacePatches() const
{
    if( !boundaryFacePatchPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelList&"
                " meshSurfaceEngine::boundaryFacePatches() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateBoundaryFacePatches();
    }

    return *boundaryFacePatchPtr_;
}

inline const labelList& meshSurfaceEngine::faceOwners() const
{
    if( !boundaryFaceOwnersPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelList& meshSurfaceEngine::faceOwners() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateBoundaryOwners();
    }

    return *boundaryFaceOwnersPtr_;
}

inline const VRWGraph& meshSurfaceEngine::pointFaces() const
{
    if( !pointFacesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::pointFaces() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculatePointFaces();
    }

    return *pointFacesPtr_;
}

inline const VRWGraph& meshSurfaceEngine::pointInFaces() const
{
    if( !pointInFacePtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::pointInFaces() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculatePointFaces();
    }

    return *pointInFacePtr_;
}
/*
inline const VRWGraph& meshSurfaceEngine::pointPatches() const
{
    if( !pointPatchesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::pointPatches() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculatePointPatches();
    }

    return *pointPatchesPtr_;
}
*/

inline const VRWGraph& meshSurfaceEngine::pointPoints() const
{
    if( !pointPointsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::pointPoints() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculatePointPoints();
    }

    return *pointPointsPtr_;
}

inline const vectorField& meshSurfaceEngine::pointNormals() const
{
    if( !pointNormalsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const vectorField& meshSurfaceEngine::pointNormals() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculatePointNormals();
    }

    return *pointNormalsPtr_;
}

inline const vectorField& meshSurfaceEngine::faceNormals() const
{
    if( !faceNormalsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const vectorField& meshSurfaceEngine::faceNormals() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateFaceNormals();
    }

    return *faceNormalsPtr_;
}

inline const vectorField& meshSurfaceEngine::faceCentres() const
{
    if( !faceCentresPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const vectorField& meshSurfaceEngine::faceCentres() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateFaceCentres();
    }

    return *faceCentresPtr_;
}

inline const edgeList& meshSurfaceEngine::edges() const
{
    if( !edgesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const edgeList& meshSurfaceEngine::edges() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateEdgesAndAddressing();
    }

    return *edgesPtr_;
}

inline const VRWGraph& meshSurfaceEngine::boundaryPointEdges() const
{
    if( !bpEdgesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::boundaryPointEdges() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateEdgesAndAddressing();
    }

    return *bpEdgesPtr_;
}

inline const VRWGraph& meshSurfaceEngine::edgeFaces() const
{
    if( !edgeFacesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::edgeFaces() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateEdgeFacesAddressing();
    }

    return *edgeFacesPtr_;
}

inline const VRWGraph& meshSurfaceEngine::faceEdges() const
{
    if( !faceEdgesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::faceEdges() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateFaceEdgesAddressing();
    }

    return *faceEdgesPtr_;
}

inline const VRWGraph& meshSurfaceEngine::edgePatches() const
{
    if( !edgePatchesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::edgePatches() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateEdgePatchesAddressing();
    }

    return *edgePatchesPtr_;
}

inline const VRWGraph& meshSurfaceEngine::faceFaces() const
{
    if( !faceFacesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::faceFaces() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateFaceFacesAddressing();
    }

    return *faceFacesPtr_;
}

inline const labelList& meshSurfaceEngine::globalBoundaryPointLabel() const
{
    if( !globalBoundaryPointLabelPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelList&"
                " meshSurfaceEngine::globalBoundaryPointLabel() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalBoundaryPointLabels();
    }

    return *globalBoundaryPointLabelPtr_;
}

inline const Map<label>&
meshSurfaceEngine::globalToLocalBndPointAddressing() const
{
    if( !globalBoundaryPointToLocalPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const Map<label>&"
                " meshSurfaceEngine::globalToLocalBndPointAddressing() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalBoundaryPointLabels();
    }

    return *globalBoundaryPointToLocalPtr_;
}

inline const VRWGraph& meshSurfaceEngine::bpAtProcs() const
{
    if( !globalBoundaryPointLabelPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::bpAtProcs() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalBoundaryPointLabels();
    }

    return *bpProcsPtr_;
}

inline const DynList<label>& meshSurfaceEngine::bpNeiProcs() const
{
    if( !bpNeiProcsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const DynList<label>& meshSurfaceEngine::bpNeiProcs() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalBoundaryPointLabels();
    }

    return *bpNeiProcsPtr_;
}

inline const labelList& meshSurfaceEngine::globalBoundaryEdgeLabel() const
{
    if( !globalBoundaryEdgeLabelPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelList&"
                " meshSurfaceEngine::globalBoundaryEdgeLabel() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalBoundaryEdgeLabels();
    }

    return *globalBoundaryEdgeLabelPtr_;
}

inline const Map<label>&
meshSurfaceEngine::globalToLocalBndEdgeAddressing() const
{
    if( !globalBoundaryEdgeToLocalPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const Map<label>&"
                " meshSurfaceEngine::globalToLocalBndEdgeAddressing() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalBoundaryEdgeLabels();
    }

    return *globalBoundaryEdgeToLocalPtr_;
}

inline const VRWGraph& meshSurfaceEngine::beAtProcs() const
{
    if( !globalBoundaryEdgeLabelPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& meshSurfaceEngine::beAtProcs() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalBoundaryEdgeLabels();
    }

    return *beProcsPtr_;
}

inline const DynList<label>& meshSurfaceEngine::beNeiProcs() const
{
    if( !beNeiProcsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const DynList<label>& meshSurfaceEngine::beNeiProcs() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalBoundaryEdgeLabels();
    }

    return *beNeiProcsPtr_;
}

inline const Map<label>& meshSurfaceEngine::otherEdgeFaceAtProc() const
{
    if( !otherEdgeFaceAtProcPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const Map<label>&"
                " meshSurfaceEngine::otherEdgeFaceAtProc() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcAddressingForProcEdges();
    }

    return *otherEdgeFaceAtProcPtr_;
}

inline const Map<label>& meshSurfaceEngine::otherEdgeFacePatch() const
{
    if( !otherEdgeFacePatchPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const Map<label>&"
                " meshSurfaceEngine::otherEdgeFacePatch() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcAddressingForProcEdges();
    }

    return *otherEdgeFacePatchPtr_;
}

inline const labelList& meshSurfaceEngine::globalBoundaryFaceLabel() const
{
    if( !globalBoundaryFaceLabelPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelList&"
                " meshSurfaceEngine::globalBoundaryFaceLabel() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalBoundaryFaceLabels();
    }

    return *globalBoundaryFaceLabelPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
