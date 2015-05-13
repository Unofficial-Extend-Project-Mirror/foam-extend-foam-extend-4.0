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

#include "cartesianMeshGenerator.H"
#include "triSurf.H"
#include "triSurfacePatchManipulator.H"
#include "demandDrivenData.H"
#include "Time.H"
#include "meshOctreeCreator.H"
#include "cartesianMeshExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceEdgeExtractorNonTopo.H"
#include "meshOptimizer.H"
#include "meshSurfaceOptimizer.H"
#include "topologicalCleaner.H"
#include "boundaryLayers.H"
#include "refineBoundaryLayers.H"
#include "renameBoundaryPatches.H"
#include "checkMeshDict.H"
#include "checkCellConnectionsOverFaces.H"
#include "checkIrregularSurfaceConnections.H"
#include "checkNonMappableCellConnections.H"
#include "checkBoundaryFacesSharingTwoEdges.H"
#include "triSurfaceMetaData.H"

//#define DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

void cartesianMeshGenerator::createCartesianMesh()
{
    //- create polyMesh from octree boxes
    cartesianMeshExtractor cme(*octreePtr_, meshDict_, mesh_);

    if( meshDict_.found("decomposePolyhedraIntoTetsAndPyrs") )
    {
        if( readBool(meshDict_.lookup("decomposePolyhedraIntoTetsAndPyrs")) )
            cme.decomposeSplitHexes();
    }

    cme.createMesh();

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_SUCCESS);
    # endif
}

void cartesianMeshGenerator::surfacePreparation()
{
    //- removes unnecessary cells and morph the boundary
    //- such that there is only one boundary face per cell
    //- It also checks topology of cells after morphing is performed
    bool changed;

    do
    {
        changed = false;

        checkIrregularSurfaceConnections checkConnections(mesh_);
        if( checkConnections.checkAndFixIrregularConnections() )
            changed = true;

        if( checkNonMappableCellConnections(mesh_).removeCells() )
            changed = true;

        if( checkCellConnectionsOverFaces(mesh_).checkCellGroups() )
            changed = true;
    } while( changed );

    checkBoundaryFacesSharingTwoEdges(mesh_).improveTopology();

    # ifdef DEBUG
    mesh_.write();
    returnReduce(1, sumOp<label>());
    ::exit(EXIT_SUCCESS);
    # endif
}

void cartesianMeshGenerator::mapMeshToSurface()
{
    //- calculate mesh surface
    meshSurfaceEngine* msePtr = new meshSurfaceEngine(mesh_);

    //- pre-map mesh surface
    meshSurfaceMapper mapper(*msePtr, *octreePtr_);
    mapper.preMapVertices();

    # ifdef DEBUG
    mesh_.write();
    returnReduce(1, sumOp<label>());
    //::exit(EXIT_SUCCESS);
    # endif

    //- map mesh surface on the geometry surface
    mapper.mapVerticesOntoSurface();

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_SUCCESS);
    # endif

    //- untangle surface faces
    meshSurfaceOptimizer(*msePtr, *octreePtr_).untangleSurface();

    # ifdef DEBUG
    mesh_.write();
    ::exit(EXIT_SUCCESS);
    # endif

    deleteDemandDrivenData(msePtr);
}

void cartesianMeshGenerator::mapEdgesAndCorners()
{
    meshSurfaceEdgeExtractorNonTopo(mesh_, *octreePtr_);

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_SUCCESS);
    # endif
}

void cartesianMeshGenerator::optimiseMeshSurface()
{
    meshSurfaceEngine mse(mesh_);
    meshSurfaceOptimizer(mse, *octreePtr_).optimizeSurface();

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_SUCCESS);
    # endif
}

void cartesianMeshGenerator::generateBoundaryLayers()
{
    //- add boundary layers
    boundaryLayers bl(mesh_);
    bl.addLayerForAllPatches();

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_SUCCESS);
    # endif
}

void cartesianMeshGenerator::refBoundaryLayers()
{
    if( meshDict_.isDict("boundaryLayers") )
    {
        refineBoundaryLayers refLayers(mesh_);

        refineBoundaryLayers::readSettings(meshDict_, refLayers);

        refLayers.refineLayers();

        meshOptimizer optimizer(mesh_);
        optimizer.untangleMeshFV();
    }
}

void cartesianMeshGenerator::optimiseFinalMesh()
{
    //- untangle the surface if needed
    bool enforceConstraints(false);
    if( meshDict_.found("enforceGeometryConstraints") )
    {
        enforceConstraints =
            readBool(meshDict_.lookup("enforceGeometryConstraints"));
    }

    if( true )
    {
        meshSurfaceEngine mse(mesh_);
        meshSurfaceOptimizer surfOpt(mse, *octreePtr_);

        if( enforceConstraints )
            surfOpt.enforceConstraints();

        surfOpt.optimizeSurface();
    }

    deleteDemandDrivenData(octreePtr_);

    //- final optimisation
    meshOptimizer optimizer(mesh_);
    if( enforceConstraints )
        optimizer.enforceConstraints();
    optimizer.optimizeMeshFV();
    optimizer.optimizeLowQualityFaces();
    optimizer.untangleMeshFV();

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_SUCCESS);
    # endif
}

void cartesianMeshGenerator::replaceBoundaries()
{
    renameBoundaryPatches rbp(mesh_, meshDict_);

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_SUCCESS);
    # endif
}

void cartesianMeshGenerator::renumberMesh()
{
    polyMeshGenModifier(mesh_).renumberMesh();

    # ifdef DEBUG
    mesh_.write();
    //::exit(EXIT_SUCCESS);
    # endif
}

void cartesianMeshGenerator::generateMesh()
{
    try
    {
        createCartesianMesh();

        surfacePreparation();

        mapMeshToSurface();

        mapEdgesAndCorners();

        optimiseMeshSurface();

        generateBoundaryLayers();

        optimiseFinalMesh();

        refBoundaryLayers();

        renumberMesh();

        replaceBoundaries();
    }
    catch(...)
    {
        WarningIn
        (
            "void cartesianMeshGenerator::generateMesh()"
        ) << "Meshing process terminated!" << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from objectRegistry
cartesianMeshGenerator::cartesianMeshGenerator(const Time& time)
:
    db_(time),
    surfacePtr_(NULL),
    meshDict_
    (
        IOobject
        (
            "meshDict",
            db_.system(),
            db_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    octreePtr_(NULL),
    mesh_(time)
{
    if( true )
    {
        checkMeshDict cmd(meshDict_);
    }

    fileName surfaceFile = meshDict_.lookup("surfaceFile");
    if( Pstream::parRun() )
        surfaceFile = ".."/surfaceFile;

    surfacePtr_ = new triSurf(db_.path()/surfaceFile);

    if( true )
    {
        //- save meta data with the mesh (surface mesh + its topology info)
        triSurfaceMetaData sMetaData(*surfacePtr_);
        const dictionary& surfMetaDict = sMetaData.metaData();

        mesh_.metaData().add("surfaceFile", surfaceFile);
        mesh_.metaData().add("surfaceMeta", surfMetaDict);
    }

    if( surfacePtr_->featureEdges().size() != 0 )
    {
        //- create surface patches based on the feature edges
        //- and update the meshDict based on the given data
        triSurfacePatchManipulator manipulator(*surfacePtr_);

        const triSurf* surfaceWithPatches =
            manipulator.surfaceWithPatches(&meshDict_);

        //- delete the old surface and assign the new one
        deleteDemandDrivenData(surfacePtr_);
        surfacePtr_ = surfaceWithPatches;
    }

    octreePtr_ = new meshOctree(*surfacePtr_);

    meshOctreeCreator(*octreePtr_, meshDict_).createOctreeBoxes();

    generateMesh();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cartesianMeshGenerator::~cartesianMeshGenerator()
{
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(octreePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void cartesianMeshGenerator::writeMesh() const
{
    mesh_.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
