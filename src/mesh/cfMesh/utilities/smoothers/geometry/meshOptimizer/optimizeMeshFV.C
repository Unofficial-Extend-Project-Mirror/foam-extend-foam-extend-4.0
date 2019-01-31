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

#include <stdexcept>

#include "demandDrivenData.H"
#include "meshOptimizer.H"
#include "polyMeshGenAddressing.H"
#include "polyMeshGenChecks.H"
#include "partTetMesh.H"
#include "HashSet.H"

#include "tetMeshOptimisation.H"
#include "boundaryLayerOptimisation.H"
#include "refineBoundaryLayers.H"
#include "meshSurfaceEngine.H"

//#define DEBUGSmooth

# ifdef DEBUGSmooth
#include "helperFunctions.H"
#include "polyMeshGenModifier.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOptimizer::untangleMeshFV
(
    const label maxNumGlobalIterations,
    const label maxNumIterations,
    const label maxNumSurfaceIterations,
    const bool relaxedCheck
)
{
    Info << "Starting untangling the mesh" << endl;

    # ifdef DEBUGSmooth
    partTetMesh tm(mesh_);
    forAll(tm.tets(), tetI)
        if( tm.tets()[tetI].mag(tm.points()) < 0.0 )
            Info << "Tet " << tetI << " is inverted!" << endl;
    polyMeshGen tetPolyMesh(mesh_.returnTime());
    tm.createPolyMesh(tetPolyMesh);
    polyMeshGenModifier(tetPolyMesh).removeUnusedVertices();
    forAll(tm.smoothVertex(), pI)
        if( !tm.smoothVertex()[pI] )
            Info << "Point " << pI << " cannot be moved!" << endl;

    const VRWGraph& pTets = tm.pointTets();
    forAll(pTets, pointI)
    {
        const LongList<partTet>& tets = tm.tets();
        forAllRow(pTets, pointI, i)
            if( tets[pTets(pointI, i)].whichPosition(pointI) < 0 )
                FatalError << "Wrong partTet" << abort(FatalError);

        partTetMeshSimplex simplex(tm, pointI);
    }

    boolList boundaryVertex(tetPolyMesh.points().size(), false);
    const labelList& neighbour = tetPolyMesh.neighbour();
    forAll(neighbour, faceI)
        if( neighbour[faceI] == -1 )
        {
            const face& f = tetPolyMesh.faces()[faceI];

            forAll(f, pI)
                boundaryVertex[f[pI]] = true;
        }

    forAll(boundaryVertex, pI)
    {
        if( boundaryVertex[pI] && tm.smoothVertex()[pI] )
            FatalErrorIn
            (
                "void meshOptimizer::untangleMeshFV()"
            ) << "Boundary vertex should not be moved!" << abort(FatalError);
    }
    # endif

    label nBadFaces, nGlobalIter(0), nIter;

    const faceListPMG& faces = mesh_.faces();

    boolList changedFace(faces.size(), true);

    //- check if any points in the tet mesh shall not move
    labelLongList lockedPoints;
    forAll(vertexLocation_, pointI)
    {
        if( vertexLocation_[pointI] & LOCKED )
            lockedPoints.append(pointI);
    }

    labelHashSet badFaces;

    do
    {
        nIter = 0;

        label minNumBadFaces(10 * faces.size()), minIter(-1);
        do
        {
            if( !relaxedCheck )
            {
                nBadFaces =
                    polyMeshGenChecks::findBadFaces
                    (
                        mesh_,
                        badFaces,
                        false,
                        &changedFace
                    );
            }
            else
            {
                nBadFaces =
                    polyMeshGenChecks::findBadFacesRelaxed
                    (
                        mesh_,
                        badFaces,
                        false,
                        &changedFace
                    );
            }

            Info << "Iteration " << nIter
                << ". Number of bad faces is " << nBadFaces << endl;

            //- perform optimisation
            if( nBadFaces == 0 )
                break;

            if( nBadFaces < minNumBadFaces )
            {
                minNumBadFaces = nBadFaces;
                minIter = nIter;
            }

            //- create a tet mesh from the mesh and the labels of bad faces
            partTetMesh tetMesh
            (
                mesh_,
                lockedPoints,
                badFaces,
                (nGlobalIter / 2) + 1
            );

            //- construct tetMeshOptimisation and improve positions of
            //- points in the tet mesh
            tetMeshOptimisation tmo(tetMesh);

            tmo.optimiseUsingKnuppMetric();

            tmo.optimiseUsingMeshUntangler();

            tmo.optimiseUsingVolumeOptimizer();

            //- update points in the mesh from the coordinates in the tet mesh
            tetMesh.updateOrigMesh(&changedFace);

        } while( (nIter < minIter+5) && (++nIter < maxNumIterations) );

        if( (nBadFaces == 0) || (++nGlobalIter >= maxNumGlobalIterations) )
            break;

        // move boundary vertices
        nIter = 0;

        while( nIter++ < maxNumSurfaceIterations )
        {
            if( !relaxedCheck )
            {
                nBadFaces =
                    polyMeshGenChecks::findBadFaces
                    (
                        mesh_,
                        badFaces,
                        false,
                        &changedFace
                    );
            }
            else
            {
                nBadFaces =
                    polyMeshGenChecks::findBadFacesRelaxed
                    (
                        mesh_,
                        badFaces,
                        false,
                        &changedFace
                    );
            }

            Info << "Iteration " << nIter
                << ". Number of bad faces is " << nBadFaces << endl;

            //- perform optimisation
            if( nBadFaces == 0 )
            {
                break;
            }
            else if( enforceConstraints_ )
            {
                const label subsetId =
                    mesh_.addPointSubset(badPointsSubsetName_);

                forAllConstIter(labelHashSet, badFaces, it)
                {
                    const face& f = faces[it.key()];
                    forAll(f, pI)
                        mesh_.addPointToSubset(subsetId, f[pI]);
                }

                WarningIn
                (
                    "void meshOptimizer::untangleMeshFV()"
                ) << "Writing mesh with " << badPointsSubsetName_
                  << " subset. These points cannot be untangled"
                  << " without sacrificing geometry constraints. Exitting.."
                  << endl;

                returnReduce(1, sumOp<label>());
                throw std::logic_error
                (
                    "void meshOptimizer::untangleMeshFV()"
                    "Cannot untangle mesh!!"
                );
            }

            //- create tethrahedral mesh from the cells which shall be smoothed
            partTetMesh tetMesh(mesh_, lockedPoints, badFaces, 0);

            //- contruct tetMeshOptimisation
            tetMeshOptimisation tmo(tetMesh);

            if( nGlobalIter < 2 )
            {
                //- the point stays in the plane determined by the point normal
                tmo.optimiseBoundaryVolumeOptimizer(true);
            }
            else if( nGlobalIter < 5 )
            {
                //- move points without any constraints on the movement
                tmo.optimiseBoundarySurfaceLaplace();
            }
            else
            {
                //- move boundary points without any constraints
                tmo.optimiseBoundaryVolumeOptimizer(false);
            }

            tetMesh.updateOrigMesh(&changedFace);

        }

    } while( nBadFaces );

    if( nBadFaces != 0 )
    {
        label subsetId = mesh_.faceSubsetIndex("badFaces");
        if( subsetId >= 0 )
            mesh_.removeFaceSubset(subsetId);
        subsetId = mesh_.addFaceSubset("badFaces");

        const labelList& owner = mesh_.owner();
        const labelList& neighbour = mesh_.neighbour();

        const label badCellsId =
                mesh_.addCellSubset("badCells");

        forAllConstIter(labelHashSet, badFaces, it)
        {
            mesh_.addFaceToSubset(subsetId, it.key());
            mesh_.addCellToSubset(badCellsId, owner[it.key()]);
            if( neighbour[it.key()] < 0 )
                continue;
            mesh_.addCellToSubset(badCellsId, neighbour[it.key()]);
        }
    }

    Info << "Finished untangling the mesh" << endl;
}

void meshOptimizer::optimizeBoundaryLayer(const bool addBufferLayer)
{
    if( mesh_.returnTime().foundObject<IOdictionary>("meshDict") )
    {
        const dictionary& meshDict =
            mesh_.returnTime().lookupObject<IOdictionary>("meshDict");

        bool smoothLayer(false);

        if( meshDict.found("boundaryLayers") )
        {
            const dictionary& layersDict = meshDict.subDict("boundaryLayers");

            if( layersDict.found("optimiseLayer") )
                smoothLayer = readBool(layersDict.lookup("optimiseLayer"));
        }

        if( !smoothLayer )
            return;

        if( addBufferLayer )
        {
            //- create a buffer layer which will not be modified by the smoother
            refineBoundaryLayers refLayers(mesh_);

            refineBoundaryLayers::readSettings(meshDict, refLayers);

            refLayers.activateSpecialMode();

            refLayers.refineLayers();

            clearSurface();
            calculatePointLocations();
        }

        Info << "Starting optimising boundary layer" << endl;

        const meshSurfaceEngine& mse = meshSurface();
        const labelList& faceOwner = mse.faceOwners();

        boundaryLayerOptimisation optimiser(mesh_, mse);

        boundaryLayerOptimisation::readSettings(meshDict, optimiser);

        optimiser.optimiseLayer();

        //- check if the bnd layer is tangled somewhere
        labelLongList bndLayerCells;
        const boolList& baseFace = optimiser.isBaseFace();

        # ifdef DEBUGSmooth
        const label blCellsId = mesh_.addCellSubset("blCells");
        # endif

        forAll(baseFace, bfI)
        {
            if( baseFace[bfI] )
            {
                bndLayerCells.append(faceOwner[bfI]);

                # ifdef DEBUGSmooth
                mesh_.addCellToSubset(blCellsId, faceOwner[bfI]);
                # endif
            }
        }

        clearSurface();
        mesh_.clearAddressingData();

        //- lock boundary layer points, faces and cells
        lockCells(bndLayerCells);

        # ifdef DEBUGSmooth
        pointField origPoints(mesh_.points().size());
        forAll(origPoints, pI)
            origPoints[pI] = mesh_.points()[pI];
        # endif

        //- optimize mesh quality
        optimizeMeshFV(5, 1, 50, 0);

        //- untangle remaining faces and lock the boundary layer cells
        untangleMeshFV(2, 50, 0);

        # ifdef DEBUGSmooth
        forAll(vertexLocation_, pI)
        {
            if( vertexLocation_[pI] & LOCKED )
            {
                if( mag(origPoints[pI] - mesh_.points()[pI]) > SMALL )
                    FatalError << "Locked points were moved"
                               << abort(FatalError);
            }
        }
        # endif

        //- unlock bnd layer points
        removeUserConstraints();

        Info << "Finished optimising boundary layer" << endl;
    }
}

void meshOptimizer::untangleBoundaryLayer()
{
    bool untangleLayer(true);
    if( mesh_.returnTime().foundObject<IOdictionary>("meshDict") )
    {
        const dictionary& meshDict =
            mesh_.returnTime().lookupObject<IOdictionary>("meshDict");

        if( meshDict.found("boundaryLayers") )
        {
            const dictionary& layersDict = meshDict.subDict("boundaryLayers");

            if( layersDict.found("untangleLayers") )
            {
                untangleLayer =
                    readBool(layersDict.lookup("untangleLayers"));
            }
        }
    }

    if( !untangleLayer )
    {
        labelHashSet badFaces;
        polyMeshGenChecks::checkFacePyramids(mesh_, false, VSMALL, &badFaces);

        const label nInvalidFaces =
            returnReduce(badFaces.size(), sumOp<label>());

        if( nInvalidFaces != 0 )
        {
            const labelList& owner = mesh_.owner();
            const labelList& neighbour = mesh_.neighbour();

            const label badBlCellsId =
                mesh_.addCellSubset("invalidBoundaryLayerCells");

            forAllConstIter(labelHashSet, badFaces, it)
            {
                mesh_.addCellToSubset(badBlCellsId, owner[it.key()]);

                if( neighbour[it.key()] < 0 )
                    continue;

                mesh_.addCellToSubset(badBlCellsId, neighbour[it.key()]);
            }

            returnReduce(1, sumOp<label>());

            throw std::logic_error
            (
                "void meshOptimizer::untangleBoundaryLayer()"
                "Found invalid faces in the boundary layer."
                " Cannot untangle mesh!!"
            );
        }
    }
    else
    {
        optimizeLowQualityFaces();
        removeUserConstraints();
        untangleMeshFV(2, 50, 1, true);
    }
}

void meshOptimizer::optimizeLowQualityFaces(const label maxNumIterations)
{
    label nBadFaces, nIter(0);

    const faceListPMG& faces = mesh_.faces();
    boolList changedFace(faces.size(), true);

    //- check if any points in the tet mesh shall not move
    labelLongList lockedPoints;
    forAll(vertexLocation_, pointI)
    {
        if( vertexLocation_[pointI] & LOCKED )
            lockedPoints.append(pointI);
    }

    do
    {
        labelHashSet lowQualityFaces;
        nBadFaces =
            polyMeshGenChecks::findLowQualityFaces
            (
                mesh_,
                lowQualityFaces,
                false,
                &changedFace
            );

        changedFace = false;
        forAllConstIter(labelHashSet, lowQualityFaces, it)
            changedFace[it.key()] = true;

        Info << "Iteration " << nIter
            << ". Number of bad faces is " << nBadFaces << endl;

        //- perform optimisation
        if( nBadFaces == 0 )
            break;

        partTetMesh tetMesh(mesh_, lockedPoints, lowQualityFaces, 2);

        //- construct tetMeshOptimisation and improve positions
        //- of points in the tet mesh
        tetMeshOptimisation tmo(tetMesh);

        tmo.optimiseUsingVolumeOptimizer();

        //- update points in the mesh from the new coordinates in the tet mesh
        tetMesh.updateOrigMesh(&changedFace);

    } while( ++nIter < maxNumIterations );
}

void meshOptimizer::optimizeMeshNearBoundaries
(
    const label maxNumIterations,
    const label numLayersOfCells
)
{
    label nIter(0);

    const faceListPMG& faces = mesh_.faces();
    boolList changedFace(faces.size(), true);

    //- check if any points in the tet mesh shall not move
    labelLongList lockedPoints;
    forAll(vertexLocation_, pointI)
    {
        if( vertexLocation_[pointI] & LOCKED )
            lockedPoints.append(pointI);
    }

    partTetMesh tetMesh(mesh_, lockedPoints, numLayersOfCells);
    tetMeshOptimisation tmo(tetMesh);
    Info << "Iteration:" << flush;
    do
    {
        tmo.optimiseUsingVolumeOptimizer(1);

        tetMesh.updateOrigMesh(&changedFace);

        Info << "." << flush;

    } while( ++nIter < maxNumIterations );

    Info << endl;
}

void meshOptimizer::optimizeMeshFV
(
    const label numLaplaceIterations,
    const label maxNumGlobalIterations,
    const label maxNumIterations,
    const label maxNumSurfaceIterations
)
{
    Info << "Starting smoothing the mesh" << endl;

    laplaceSmoother lps(mesh_, vertexLocation_);
    lps.optimizeLaplacianPC(numLaplaceIterations);

    untangleMeshFV
    (
        maxNumGlobalIterations,
        maxNumIterations,
        maxNumSurfaceIterations
    );

    Info << "Finished smoothing the mesh" << endl;
}

void meshOptimizer::optimizeMeshFVBestQuality
(
    const label maxNumIterations,
    const scalar threshold
)
{
    label nBadFaces, nIter(0);
    label minIter(-1);

    const faceListPMG& faces = mesh_.faces();
    boolList changedFace(faces.size(), true);

    //- check if any points in the tet mesh shall not move
    labelLongList lockedPoints;
    forAll(vertexLocation_, pointI)
    {
        if( vertexLocation_[pointI] & LOCKED )
            lockedPoints.append(pointI);
    }

    label minNumBadFaces(10 * faces.size());
    do
    {
        labelHashSet lowQualityFaces;
        nBadFaces =
            polyMeshGenChecks::findWorstQualityFaces
            (
                mesh_,
                lowQualityFaces,
                false,
                &changedFace,
                threshold
            );

        changedFace = false;
        forAllConstIter(labelHashSet, lowQualityFaces, it)
            changedFace[it.key()] = true;

        Info << "Iteration " << nIter
            << ". Number of worst quality faces is " << nBadFaces << endl;

        //- perform optimisation
        if( nBadFaces == 0 )
            break;

        if( nBadFaces < minNumBadFaces )
        {
            minNumBadFaces = nBadFaces;

            //- update the iteration number when the minimum is achieved
            minIter = nIter;
        }

        partTetMesh tetMesh(mesh_, lockedPoints, lowQualityFaces, 2);

        //- construct tetMeshOptimisation and improve positions
        //- of points in the tet mesh
        tetMeshOptimisation tmo(tetMesh);

        tmo.optimiseUsingVolumeOptimizer(20);

        //- update points in the mesh from the new coordinates in the tet mesh
        tetMesh.updateOrigMesh(&changedFace);

    } while( (nIter < minIter+5) && (++nIter < maxNumIterations) );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
