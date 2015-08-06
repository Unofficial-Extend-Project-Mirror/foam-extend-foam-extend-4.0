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
#include "meshSurfaceEngineModifier.H"
#include "meshSurfaceMapper2D.H"
#include "polyMeshGen2DEngine.H"
#include "meshOctree.H"
#include "refLabelledPoint.H"
#include "helperFunctionsPar.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper2D::preMapVertices(const label nIterations)
{
    Info << "Smoothing mesh surface before mapping. Iteration:" << flush;

    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bp = surfaceEngine_.bp();
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const VRWGraph& edgeFaces = surfaceEngine_.edgeFaces();
    const edgeList& edges = surfaceEngine_.edges();

    List<labelledPoint> preMapPositions(activeBoundaryEdges_.size());

    for(label iterI=0;iterI<nIterations;++iterI)
    {
        labelLongList parBndEdges;

        //- use the shrinking laplace first
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 40)
        # endif
        forAll(activeBoundaryEdges_, eI)
        {
            const label beI = activeBoundaryEdges_[eI];

            labelledPoint lp(0, vector::zero);

            if( edgeFaces.sizeOfRow(beI) == 2 )
            {
                forAllRow(edgeFaces, beI, efI)
                {
                    ++lp.pointLabel();
                    lp.coordinates() += faceCentres[edgeFaces(beI, efI)];
                }

            }
            else if( edgeFaces.sizeOfRow(beI) == 1 )
            {
                ++lp.pointLabel();
                lp.coordinates() += faceCentres[edgeFaces(beI, 0)];

                #ifdef USE_OMP
                # pragma omp critical
                # endif
                parBndEdges.append(eI);
            }

            //- store the information
            preMapPositions[eI] = lp;
        }

        if( Pstream::parRun() )
        {
            const VRWGraph& beAtProcs = surfaceEngine_.beAtProcs();
            const labelList& globalEdgeLabel =
                surfaceEngine_.globalBoundaryEdgeLabel();
            const Map<label>& globalToLocal =
                surfaceEngine_.globalToLocalBndEdgeAddressing();

            //- collect data to be sent to other processors
            std::map<label, LongList<refLabelledPoint> > exchangeData;
            forAll(surfaceEngine_.beNeiProcs(), i)
                exchangeData.insert
                (
                    std::make_pair
                    (
                        surfaceEngine_.beNeiProcs()[i],
                        LongList<refLabelledPoint>()
                    )
                );

            Map<label> edgeToActiveAddressing;
            forAll(parBndEdges, i)
            {
                const label beI = activeBoundaryEdges_[parBndEdges[i]];
                edgeToActiveAddressing.insert(beI, parBndEdges[i]);

                forAllRow(beAtProcs, beI, procI)
                {
                    const label neiProc = beAtProcs(beI, procI);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append
                    (
                        refLabelledPoint
                        (
                            globalEdgeLabel[beI],
                            preMapPositions[beI]
                        )
                    );
                }
            }

            //- exchange data with other processors
            LongList<refLabelledPoint> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            //- combine collected data with the available data
            forAll(receivedData, i)
            {
                const refLabelledPoint& rlp = receivedData[i];
                const labelledPoint& lps = rlp.lPoint();

                const label beI = globalToLocal[rlp.objectLabel()];
                const label eI = edgeToActiveAddressing[beI];

                labelledPoint& lp = preMapPositions[eI];
                lp.pointLabel() += lps.pointLabel();
                lp.coordinates() += lps.coordinates();
            }
        }

        //- calculate coordinates of points for searching
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(activeBoundaryEdges_, eI)
        {
            labelledPoint& lp = preMapPositions[eI];

            if( lp.pointLabel() == 0 )
            {
                Warning << "Surface edge " << activeBoundaryEdges_[eI]
                    << " has no active faces" << endl;
                continue;
            }

            lp.coordinates() /= lp.pointLabel();
        }

        //- create the surface modifier and move the surface points
        meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(activeBoundaryEdges_, eI)
        {
            const label beI = activeBoundaryEdges_[eI];
            const edge& e = edges[beI];

            const point& p = points[e.start()];

            label patch, nt;
            point pMap = p;
            scalar dSq;

            meshOctree_.findNearestSurfacePoint
            (
                pMap,
                dSq,
                nt,
                patch,
                preMapPositions[eI].coordinates()
            );

            pMap.z() = p.z();
            point newP = 0.5 * (pMap + p);

            surfaceModifier.moveBoundaryVertexNoUpdate(bp[e.start()], newP);
            newP.z() = points[e.end()].z();
            surfaceModifier.moveBoundaryVertexNoUpdate(bp[e.end()], newP);
        }

        surfaceModifier.updateGeometry();
        surfaceModifier.syncVerticesAtParallelBoundaries();

        Info << "." << flush;
    }

    Info << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
