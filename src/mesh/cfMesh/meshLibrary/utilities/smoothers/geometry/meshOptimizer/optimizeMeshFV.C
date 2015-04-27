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
#include "partTetMesh.H"
#include "HashSet.H"

#include "tetMeshOptimisation.H"

//#define DEBUGSmooth

# ifdef DEBUGSmooth
#include "helperFunctions.H"
#include "polyMeshGenModifier.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOptimizer::untangleMeshFV()
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
    const label maxNumGlobalIterations(10);

    const faceListPMG& faces = mesh_.faces();
    boolList changedFace(faces.size(), true);

    labelHashSet badFaces;

    do
    {
        nIter = 0;

        label minNumBadFaces(10 * faces.size()), minIter(-1);
        do
        {
            nBadFaces = findBadFaces(badFaces, changedFace);

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

            partTetMesh tetMesh(mesh_, badFaces, (nGlobalIter / 5) + 1);

            tetMeshOptimisation tmo(tetMesh);

            tmo.optimiseUsingKnuppMetric();

            tmo.optimiseUsingMeshUntangler();

            tmo.optimiseUsingVolumeOptimizer();

            tetMesh.updateOrigMesh(&changedFace);

        } while( (nIter < minIter+5) && (++nIter < 50) );

        if( (nBadFaces == 0) || (++nGlobalIter >= maxNumGlobalIterations) )
            break;

        // move boundary vertices
        nIter = 0;

        do
        {
            nBadFaces = findBadFaces(badFaces, changedFace);

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

            partTetMesh tetMesh(mesh_, badFaces, 0);
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

        } while( ++nIter < 2 );
    }
    while( nBadFaces );

    Info << "Finished untangling the mesh" << endl;
}

void meshOptimizer::optimizeLowQualityFaces()
{
    label nBadFaces, nIter(0);

    const faceListPMG& faces = mesh_.faces();
    boolList changedFace(faces.size(), true);

    label minNumBadFaces(10 * faces.size()), minIter(-1);
    do
    {
        labelHashSet lowQualityFaces;
        nBadFaces = findLowQualityFaces(lowQualityFaces, changedFace);

        changedFace = false;
        forAllConstIter(labelHashSet, lowQualityFaces, it)
        changedFace[it.key()] = true;

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

        partTetMesh tetMesh(mesh_, lowQualityFaces, 2);

        tetMeshOptimisation tmo(tetMesh);

        tmo.optimiseUsingKnuppMetric();

        tmo.optimiseUsingMeshUntangler();

        tmo.optimiseUsingVolumeOptimizer();

        tetMesh.updateOrigMesh(&changedFace);

    } while( (nIter < minIter+2) && (++nIter < 10) );
}

void meshOptimizer::optimizeMeshFV()
{
    Info << "Starting smoothing the mesh" << endl;

    laplaceSmoother lps(mesh_, vertexLocation_);
    lps.optimizeLaplacianPC(5);

    untangleMeshFV();

    Info << "Finished smoothing the mesh" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
