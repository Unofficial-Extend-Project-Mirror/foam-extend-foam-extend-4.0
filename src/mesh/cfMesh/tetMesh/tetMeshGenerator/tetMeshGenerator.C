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

#include "tetMeshGenerator.H"
#include "triSurf.H"
#include "demandDrivenData.H"
#include "Time.H"
#include "meshOctreeCreator.H"
#include "meshOctreeAutomaticRefinement.H"
#include "tetMeshExtractorOctree.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceEdgeExtractorNonTopo.H"
#include "surfaceMorpherCells.H"
#include "meshOptimizer.H"
#include "meshSurfaceOptimizer.H"
#include "topologicalCleaner.H"
#include "boundaryLayers.H"
#include "renameBoundaryPatches.H"
#include "checkMeshDict.H"
#include "triSurfacePatchManipulator.H"
#include "refineBoundaryLayers.H"
#include "triSurfaceMetaData.H"

//#define DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

void tetMeshGenerator::createTetMesh()
{
    //- create tet Mesh from octree and Delaunay tets
    tetMeshExtractorOctree tme(*octreePtr_, meshDict_, mesh_);

    tme.createMesh();

    # ifdef DEBUG
    mesh_.write();
    //::exit(0);
    # endif
}

void tetMeshGenerator::surfacePreparation()
{
    //- removes unnecessary cells and morph the boundary
    //- such that there is only one boundary face per cell
    //- It also checks topology of cells after morphing is performed
    do
    {
        surfaceMorpherCells* cmPtr = new surfaceMorpherCells(mesh_);
        cmPtr->morphMesh();
        deleteDemandDrivenData(cmPtr);
    }
    while( topologicalCleaner(mesh_).cleanTopology() );

    # ifdef DEBUG
    mesh_.write();
    //::exit(0);
    # endif
}

void tetMeshGenerator::mapMeshToSurface()
{
    //- calculate mesh surface
    meshSurfaceEngine* msePtr = new meshSurfaceEngine(mesh_);

    //- map mesh surface on the geometry surface
    meshSurfaceMapper(*msePtr, *octreePtr_).mapVerticesOntoSurface();

    # ifdef DEBUG
    mesh_.write();
    //::exit(0);
    # endif

    //- untangle surface faces
    meshSurfaceOptimizer(*msePtr, *octreePtr_).untangleSurface();

    # ifdef DEBUG
    mesh_.write();
    //::exit(0);
    # endif

    deleteDemandDrivenData(msePtr);
}

void tetMeshGenerator::mapEdgesAndCorners()
{
    meshSurfaceEdgeExtractorNonTopo(mesh_, *octreePtr_);

    # ifdef DEBUG
    mesh_.write();
    //::exit(0);
    # endif
}

void tetMeshGenerator::optimiseMeshSurface()
{
    meshSurfaceEngine mse(mesh_);
    meshSurfaceOptimizer(mse, *octreePtr_).optimizeSurface();

    # ifdef DEBUG
    mesh_.write();
    //::exit(0);
    # endif
}

void tetMeshGenerator::generateBoudaryLayers()
{
    if( meshDict_.found("boundaryLayers") )
    {
        boundaryLayers bl(mesh_);

        const dictionary& bndLayers = meshDict_.subDict("boundaryLayers");

        if( bndLayers.found("nLayers") )
        {
            const label nLayers = readLabel(bndLayers.lookup("nLayers"));

            if( nLayers > 0 )
                bl.addLayerForAllPatches();
        }
        else if( bndLayers.found("patchBoundaryLayers") )
        {
            const dictionary& patchLayers =
                bndLayers.subDict("patchBoundaryLayers");
            const wordList createLayers = patchLayers.toc();

            forAll(createLayers, patchI)
                bl.addLayerForPatch(createLayers[patchI]);
        }
    }

    # ifdef DEBUG
    mesh_.write();
    //::exit(0);
    # endif
}

void tetMeshGenerator::optimiseFinalMesh()
{
    //- final optimisation
    bool enforceConstraints(false);
    if( meshDict_.found("enforceGeometryConstraints") )
    {
        enforceConstraints =
            readBool(meshDict_.lookup("enforceGeometryConstraints"));
    }

    meshOptimizer optimizer(mesh_);
    if( enforceConstraints )
        optimizer.enforceConstraints();

    optimizer.optimizeSurface(*octreePtr_);

    deleteDemandDrivenData(octreePtr_);

    optimizer.optimizeMeshFV();
    optimizer.optimizeLowQualityFaces();
    optimizer.optimizeMeshFV();

    # ifdef DEBUG
    mesh_.write();
    //::exit(0);
    # endif
}

void tetMeshGenerator::refBoundaryLayers()
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

void tetMeshGenerator::replaceBoundaries()
{
    renameBoundaryPatches rbp(mesh_, meshDict_);

    # ifdef DEBUG
    mesh_.write();
    //::exit(0);
    # endif
}

void tetMeshGenerator::renumberMesh()
{
    polyMeshGenModifier(mesh_).renumberMesh();

    # ifdef DEBUG
    mesh_.write();
    //::exit(0);
    # endif
}

void tetMeshGenerator::generateMesh()
{
    try
    {
        createTetMesh();

        surfacePreparation();

        mapMeshToSurface();

        mapEdgesAndCorners();

        optimiseMeshSurface();

        generateBoudaryLayers();

        optimiseFinalMesh();

        refBoundaryLayers();

        renumberMesh();

        replaceBoundaries();
    }
    catch(...)
    {
        WarningIn
        (
            "void tetMeshGenerator::generateMesh()"
        ) << "Meshing process terminated!" << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from Time
tetMeshGenerator::tetMeshGenerator(const Time& time)
:
    runTime_(time),
    surfacePtr_(NULL),
    meshDict_
    (
        IOobject
        (
            "meshDict",
            runTime_.system(),
            runTime_,
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

    const fileName surfaceFile = meshDict_.lookup("surfaceFile");

    surfacePtr_ = new triSurf(runTime_.path()/surfaceFile);

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

    meshOctreeCreator* octreeCreatorPtr =
        new meshOctreeCreator(*octreePtr_, meshDict_);
    octreeCreatorPtr->createOctreeBoxes();
    deleteDemandDrivenData(octreeCreatorPtr);

    generateMesh();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tetMeshGenerator::~tetMeshGenerator()
{
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(octreePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetMeshGenerator::writeMesh() const
{
    mesh_.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
