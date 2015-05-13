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
#include "tetMeshOptimisation.H"
#include "partTetMesh.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"

#include "partTetMeshSimplex.H"
#include "meshUntangler.H"
#include "volumeOptimizer.H"
#include "knuppMetric.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
tetMeshOptimisation::tetMeshOptimisation(partTetMesh& mesh)
:
    tetMesh_(mesh)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tetMeshOptimisation::~tetMeshOptimisation()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetMeshOptimisation::optimiseUsingKnuppMetric()
{
    const LongList<point>& points = tetMesh_.points();
    const LongList<partTet>& tets = tetMesh_.tets();
    const LongList<direction>& smoothVertex = tetMesh_.smoothVertex();

    boolList negativeNode(smoothVertex.size()), invertedTets(tets.size());

    //- try getting rid of negative volume using the Patrik Knupp's metric
    //- which gets non-negative contributions from invertex tets, only
    # ifdef USE_OMP
    # pragma omp parallel for if( tets.size() > 100 ) \
    schedule(dynamic, 10)
    # endif
    forAll(tets, tetI)
    {
        invertedTets[tetI] = false;

        if( tets[tetI].mag(points) < VSMALL )
            invertedTets[tetI] = true;
    }

    label nIter(0), nNegative, nNegativeBefore;

    do
    {
        //- find the number of inverted tets
        nNegative = 0;
        negativeNode = false;
        # ifdef USE_OMP
        # pragma omp parallel for if( tets.size() > 100 ) \
        schedule(dynamic, 10) reduction(+ : nNegative)
        # endif
        forAll(invertedTets, tetI)
        {
            if( invertedTets[tetI] )
            {
                ++nNegative;
                const partTet& tet = tets[tetI];

                for(label i=0;i<4;++i)
                    negativeNode[tet[i]] = true;
            }
        }

        reduce(nNegative, sumOp<label>());
        if( nNegative == 0 )
            return;

        //- make sure that the points at procesor boundaries are selected
        //- at all processors
        if( Pstream::parRun() )
            unifyNegativePoints(negativeNode);

        //- smooth the mesh
        List<LongList<labelledPoint> > newPositions;
        # ifdef USE_OMP
        # pragma omp parallel if( smoothVertex.size() > 100 )
        # endif
        {
            # ifdef USE_OMP
            # pragma omp master
            {
                newPositions.setSize(omp_get_num_threads());
            }

            # pragma omp barrier

            LongList<labelledPoint>& np = newPositions[omp_get_thread_num()];
            # else
            newPositions.setSize(1);
            LongList<labelledPoint>& np = newPositions[0];
            # endif

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 10)
            # endif
            forAll(smoothVertex, nodeI)
            {
                if( !negativeNode[nodeI] )
                    continue;

                if( smoothVertex[nodeI] & partTetMesh::SMOOTH )
                {
                    partTetMeshSimplex simplex(tetMesh_, nodeI);
                    knuppMetric(simplex).optimizeNodePosition();
                    np.append(labelledPoint(nodeI, simplex.centrePoint()));
                }
            }
        }

        //- update mesh vertices
        tetMesh_.updateVerticesSMP(newPositions);
        newPositions.clear();

        if( Pstream::parRun() )
        {
            updateBufferLayerPoints();
            unifyCoordinatesParallel(&negativeNode);
        }

        //- check which tets have been repaired
        boolList helper(invertedTets.size());
        nNegativeBefore = nNegative;
        nNegative = 0;

        # ifdef USE_OMP
        # pragma omp parallel for if( tets.size() > 100 ) \
        schedule(dynamic, 10) reduction(+ : nNegative)
        # endif
        forAll(tets, tetI)
        {
            helper[tetI] = false;

            if( invertedTets[tetI] && (tets[tetI].mag(points) < VSMALL) )
            {
                helper[tetI] = true;

                const partTet& tet = tets[tetI];

                for(label i=0;i<4;++i)
                    negativeNode[tet[i]] = true;
            }
        }
        invertedTets.transfer(helper);

        reduce(nNegative, sumOp<label>());
        if( nNegative == 0 )
            return;

    } while( (nNegative < nNegativeBefore) || (++nIter < 5) );
}

void tetMeshOptimisation::optimiseUsingMeshUntangler()
{
    const LongList<point>& points = tetMesh_.points();
    const LongList<partTet>& tets = tetMesh_.tets();
    const LongList<direction>& smoothVertex = tetMesh_.smoothVertex();

    boolList negativeNode(smoothVertex.size()), invertedTets(tets.size());

    //- try getting rid of negative volume using the untangler
    # ifdef USE_OMP
    # pragma omp parallel for if( tets.size() > 100 ) \
    schedule(dynamic, 10)
    # endif
    forAll(tets, tetI)
    {
        invertedTets[tetI] = false;

        if( tets[tetI].mag(points) < VSMALL )
            invertedTets[tetI] = true;
    }

    label nIter(0), nNegative, nNegativeBefore;

    do
    {
        //- find the number of inverted tets
        nNegative = 0;
        negativeNode = false;
        # ifdef USE_OMP
        # pragma omp parallel for if( tets.size() > 100 ) \
        schedule(dynamic, 10) reduction(+ : nNegative)
        # endif
        forAll(invertedTets, tetI)
        {
            if( invertedTets[tetI] )
            {
                ++nNegative;
                const partTet& tet = tets[tetI];

                for(label i=0;i<4;++i)
                    negativeNode[tet[i]] = true;
            }
        }

        reduce(nNegative, sumOp<label>());
        if( nNegative == 0 )
            return;

        //- make sure that the points at procesor boundaries are selected
        //- at all processors
        if( Pstream::parRun() )
            unifyNegativePoints(negativeNode);

        //- smooth the mesh
        List<LongList<labelledPoint> > newPositions;
        # ifdef USE_OMP
        # pragma omp parallel if( smoothVertex.size() > 100 )
        # endif
        {
            # ifdef USE_OMP
            # pragma omp master
            {
                newPositions.setSize(omp_get_num_threads());
            }

            # pragma omp barrier

            LongList<labelledPoint>& np = newPositions[omp_get_thread_num()];
            # else
            newPositions.setSize(1);
            LongList<labelledPoint>& np = newPositions[0];
            # endif

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 10)
            # endif
            forAll(smoothVertex, nodeI)
            {
                if( !negativeNode[nodeI] )
                    continue;

                if( smoothVertex[nodeI] & partTetMesh::SMOOTH )
                {
                    partTetMeshSimplex simplex(tetMesh_, nodeI);
                    meshUntangler(simplex).optimizeNodePosition();
                    np.append(labelledPoint(nodeI, simplex.centrePoint()));
                }
            }
        }

        //- update mesh vertices
        tetMesh_.updateVerticesSMP(newPositions);
        newPositions.clear();

        if( Pstream::parRun() )
        {
            updateBufferLayerPoints();
            unifyCoordinatesParallel(&negativeNode);
        }

        //- check which tets have been repaired
        boolList helper(invertedTets.size());
        nNegativeBefore = nNegative;
        nNegative = 0;
        # ifdef USE_OMP
        # pragma omp parallel for if( tets.size() > 100 ) \
        schedule(dynamic, 10) reduction(+ : nNegative)
        # endif
        forAll(tets, tetI)
        {
            helper[tetI] = false;

            if( invertedTets[tetI] && (tets[tetI].mag(points) < VSMALL) )
            {
                ++nNegative;
                const partTet& tet = tets[tetI];

                for(label i=0;i<4;++i)
                    negativeNode[tet[i]] = true;
            }
        }

        reduce(nNegative, sumOp<label>());
        if( nNegative == 0 )
            return;
        invertedTets.transfer(helper);

    } while( (nNegative < nNegativeBefore) || (++nIter < 5) );
}

void tetMeshOptimisation::optimiseUsingVolumeOptimizer()
{
    const LongList<direction>& smoothVertex = tetMesh_.smoothVertex();

    //- use mesh optimizer to improve the result
    for(label i=0;i<2;++i)
    {
        List<LongList<labelledPoint> > newPositions;

        # ifdef USE_OMP
        # pragma omp parallel if( smoothVertex.size() > 100 )
        # endif
        {
            # ifdef USE_OMP
            # pragma omp master
            {
                newPositions.setSize(omp_get_num_threads());
            }

            # pragma omp barrier

            LongList<labelledPoint>& np = newPositions[omp_get_thread_num()];
            # else
            newPositions.setSize(1);
            LongList<labelledPoint>& np = newPositions[0];
            # endif

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 10)
            # endif
            forAll(smoothVertex, nodeI)
                if( smoothVertex[nodeI] & partTetMesh::SMOOTH )
                {
                    partTetMeshSimplex simplex(tetMesh_, nodeI);

                    volumeOptimizer vOpt(simplex);
                    vOpt.optimizeNodePosition(1e-5);

                    np.append(labelledPoint(nodeI, simplex.centrePoint()));
                }
        }

        //- update mesh vertices
        tetMesh_.updateVerticesSMP(newPositions);
        newPositions.clear();

        if( Pstream::parRun() )
        {
            updateBufferLayerPoints();
            unifyCoordinatesParallel();
        }
    }
}

void tetMeshOptimisation::optimiseBoundaryVolumeOptimizer
(
    const bool nonShrinking
)
{
    const LongList<point>& points = tetMesh_.points();
    const LongList<direction>& smoothVertex = tetMesh_.smoothVertex();

    # ifdef USE_OMP
    label nThreads = omp_get_num_procs();
    if( smoothVertex.size() < 100 )
      nThreads = 1;
    # else
    const label nThreads(1);
    # endif

    for(label i=0;i<3;++i)
    {
        List<LongList<labelledPoint> > newPositions(nThreads);

        # ifdef USE_OMP
        # pragma omp parallel num_threads(nThreads)
        # endif
        {
            # ifdef USE_OMP
            LongList<labelledPoint>& np = newPositions[omp_get_thread_num()];
            # else
            LongList<labelledPoint>& np = newPositions[0];
            # endif

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 5)
            # endif
            forAll(smoothVertex, nodeI)
                if( smoothVertex[nodeI] & partTetMesh::BOUNDARY )
                {
                    partTetMeshSimplex simplex(tetMesh_, nodeI);

                    volumeOptimizer vOpt(simplex);
                    vOpt.optimizeNodePosition(1e-5);

                    if( nonShrinking )
                    {
                        //- find boundary faces of the simplex
                        const DynList<point, 128>& pts = simplex.pts();
                        const DynList<partTet, 128>& tets = simplex.tets();
                        DynList<edge, 64> bEdges;
                        DynList<label, 64> numAppearances;

                        forAll(tets, tetI)
                        {
                            const partTet& tet = tets[tetI];
                            for(label i=0;i<3;++i)
                            {
                                edge e(tet[i], tet[(i+1)%3]);

                                const label pos = bEdges.containsAtPosition(e);
                                if( pos < 0 )
                                {
                                    bEdges.append(e);
                                    numAppearances.append(1);
                                }
                                else
                                {
                                    ++numAppearances(pos);
                                }
                            }
                        }

                        //- create normal tensor of the simplex
                        symmTensor nt(symmTensor::zero);
                        forAll(bEdges, beI)
                        {
                            if( numAppearances[beI] != 1 )
                                continue;

                            triangle<point, point> tri
                            (
                                pts[bEdges[beI].start()],
                                pts[bEdges[beI].end()],
                                points[nodeI]
                            );

                            vector n = tri.normal();
                            n /= (mag(n) + VSMALL);
                            for(direction i=0;i<vector::nComponents;++i)
                                if( Foam::mag(n[i]) < (100. * SMALL) )
                                    n[i] = 0.0;

                            nt += symm(n * n);
                        }

                        const vector ev = eigenValues(nt);

                        //- make sure the point stays on the surface
                        vector disp = simplex.centrePoint() - points[nodeI];

                        if( mag(ev[2]) > (mag(ev[1]) + mag(ev[0])) )
                        {
                            //- ordinary surface vertex
                            vector normal = eigenVector(nt, ev[2]);
                            normal /= (mag(normal)+VSMALL);
                            disp -= (disp & normal) * normal;
                        }
                        else if( mag(ev[1]) > 0.5 * (mag(ev[2]) + mag(ev[0])) )
                        {
                            //- this vertex is on an edge
                            vector normal1 = eigenVector(nt, ev[1]);
                            normal1 /= (mag(normal1)+VSMALL);
                            vector normal2 = eigenVector(nt, ev[2]);
                            normal2 /= (mag(normal2)+VSMALL);

                            vector eVec = normal1 ^ normal2;
                            eVec /= (mag(eVec) + VSMALL);

                            disp = (disp & eVec) * eVec;
                        }
                        else
                        {
                            //- this vertex is a corner. do not move it
                            continue;
                        }

                        const point newP = points[nodeI] + disp;
                        np.append(labelledPoint(nodeI, newP));
                    }
                    else
                    {
                        //- move the vertex without constraining it
                        np.append(labelledPoint(nodeI, simplex.centrePoint()));
                    }
                }
        }

        //- update tetMesh
        tetMesh_.updateVerticesSMP(newPositions);
        newPositions.clear();

        if( Pstream::parRun() )
        {
            updateBufferLayerPoints();
            unifyCoordinatesParallel();
        }
    }
}

void tetMeshOptimisation::optimiseBoundarySurfaceLaplace()
{
    const LongList<direction>& smoothVertex = tetMesh_.smoothVertex();

    # ifdef USE_OMP
    label nThreads = omp_get_num_procs();
    if( smoothVertex.size() < 1000 )
      nThreads = 1;
    # else
    const label nThreads(1);
    # endif

    for(label i=0;i<3;++i)
    {
        List<LongList<labelledPoint> > newPositions(nThreads);

        # ifdef USE_OMP
        # pragma omp parallel num_threads(nThreads)
        # endif
        {
            # ifdef USE_OMP
            LongList<labelledPoint>& np =
                newPositions[omp_get_thread_num()];
            # else
            LongList<labelledPoint>& np = newPositions[0];
            # endif

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 5)
            # endif
            forAll(smoothVertex, nodeI)
            {
                if( smoothVertex[nodeI] & partTetMesh::BOUNDARY )
                {
                    partTetMeshSimplex simplex(tetMesh_, nodeI);

                    //- find boundary faces of the simplex
                    const DynList<point, 128>& pts = simplex.pts();
                    const DynList<partTet, 128>& tets = simplex.tets();
                    DynList<edge, 64> bndEdges;
                    DynList<label, 64> numAppearances;

                    //- find boundary edges of the simplex
                    forAll(tets, tetI)
                    {
                        const partTet& tet = tets[tetI];
                        for(label i=0;i<3;++i)
                        {
                            const edge e(tet[i], tet[(i+1)%3]);
                            const label pos = bndEdges.containsAtPosition(e);

                            if( pos < 0 )
                            {
                                bndEdges.append(e);
                                numAppearances.append(1);
                            }
                            else
                            {
                                ++numAppearances(pos);
                            }
                        }
                    }

                    point newP(vector::zero);
                    label counter(0);
                    forAll(bndEdges, beI)
                    {
                        if( numAppearances[beI] != 1 )
                            continue;

                        triangle<point, point> tri
                        (
                            pts[bndEdges[beI].start()],
                            pts[bndEdges[beI].end()],
                            simplex.centrePoint()
                        );

                        newP += tri.centre();
                        ++counter;
                    }

                    if( counter != 0 )
                    {
                        newP /= counter;
                        np.append(labelledPoint(nodeI, newP));
                    }
                }
            }
        }

        //- update tetMesh with new vertex positions
        tetMesh_.updateVerticesSMP(newPositions);
        newPositions.clear();

        if( Pstream::parRun() )
        {
            updateBufferLayerPoints();
            unifyCoordinatesParallel();
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
