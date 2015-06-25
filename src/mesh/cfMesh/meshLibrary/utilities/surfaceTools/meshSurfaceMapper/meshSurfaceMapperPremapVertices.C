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
#include "meshSurfaceMapper.H"
#include "meshOctree.H"
#include "triSurf.H"
#include "refLabelledPoint.H"
#include "refLabelledPointScalar.H"
#include "helperFunctions.H"
#include "meshSurfaceOptimizer.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper::preMapVertices(const label nIterations)
{
    Info << "Smoothing mesh surface before mapping." << endl;

    const labelList& boundaryPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();
    const VRWGraph& pointInFace = surfaceEngine_.pointInFaces();
    const labelList& bp = surfaceEngine_.bp();
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

    const triSurf& surf = meshOctree_.surface();

    List<labelledPointScalar> preMapPositions(boundaryPoints.size());
    List<DynList<scalar, 6> > faceCentreDistances(bFaces.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(bFaces, bfI)
    {
        const point& c = faceCentres[bfI];
        const face& bf = bFaces[bfI];

        faceCentreDistances[bfI].setSize(bf.size());

        forAll(bf, pI)
        {
            faceCentreDistances[bfI][pI] = magSqr(points[bf[pI]] - c);
        }
    }

    for(label iterI=0;iterI<nIterations;++iterI)
    {
        //- find patches in the vicinity of a boundary face
        List<DynList<label> > boundaryPointPatches(boundaryPoints.size());
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];
            scalar boxSize(0.0);
            forAll(bf, pI)
            {
                boxSize =
                    Foam::max
                    (
                        boxSize,
                        mag(faceCentres[bfI] - points[bf[pI]])
                    );
            }

            const boundBox bb
            (
                faceCentres[bfI] - vector(boxSize, boxSize, boxSize),
                faceCentres[bfI] + vector(boxSize, boxSize, boxSize)
            );

            DynList<label> containedLeaves;
            meshOctree_.findLeavesContainedInBox(bb, containedLeaves);

            DynList<label> patches;
            forAll(containedLeaves, clI)
            {
                DynList<label> ct;
                meshOctree_.containedTriangles(containedLeaves[clI], ct);

                forAll(ct, i)
                    patches.appendIfNotIn(surf[ct[i]].region());
            }

            scalar metric(VGREAT);
            label bestPatch(-1);
            forAll(patches, ptchI)
            {
                const scalar m = faceMetricInPatch(bfI, patches[ptchI]);

                if( m < metric )
                {
                    metric = m;
                    bestPatch = patches[ptchI];
                }
            }

            forAll(bf, pI)
                boundaryPointPatches[bp[bf[pI]]].appendIfNotIn(bestPatch);
        }

        //- use the shrinking laplace first
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 40)
        # endif
        forAll(pointFaces, bpI)
        {
            labelledPointScalar lp(bpI, vector::zero, 0.0);

            const point& p = points[boundaryPoints[bpI]];

            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);
                const point& fc = faceCentres[pointFaces(bpI, pfI)];
                const label pos = pointInFace(bpI, pfI);
                const scalar w
                (
                    max(magSqr(p - fc) / faceCentreDistances[bfI][pos], SMALL)
                );
                lp.coordinates() += w * faceCentres[bfI];
                lp.scalarValue() += w;
            }

            preMapPositions[bpI] = lp;
        }

        //- pointer needed in case of parallel calculation
        const VRWGraph* bpAtProcsPtr(NULL);

        if( Pstream::parRun() )
        {
            const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
            bpAtProcsPtr = &bpAtProcs;
            const labelList& globalPointLabel =
                surfaceEngine_.globalBoundaryPointLabel();
            const Map<label>& globalToLocal =
                surfaceEngine_.globalToLocalBndPointAddressing();

            //- collect data to be sent to other processors
            std::map<label, LongList<labelledPointScalar> > exchangeData;
            forAll(surfaceEngine_.bpNeiProcs(), i)
                exchangeData.insert
                (
                    std::make_pair
                    (
                        surfaceEngine_.bpNeiProcs()[i],
                        LongList<labelledPointScalar>()
                    )
                );

            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                forAllRow(bpAtProcs, bpI, procI)
                {
                    const label neiProc = bpAtProcs(bpI, procI);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append
                    (
                        labelledPointScalar
                        (
                            globalPointLabel[bpI],
                            preMapPositions[bpI].coordinates(),
                            preMapPositions[bpI].scalarValue()
                        )
                    );
                }
            }

            //- exchange data with other processors
            LongList<labelledPointScalar> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            //- combine collected data with the available data
            forAll(receivedData, i)
            {
                const labelledPointScalar& lps = receivedData[i];

                const label bpI = globalToLocal[lps.pointLabel()];

                labelledPointScalar& lp = preMapPositions[bpI];
                lp.coordinates() += lps.coordinates();
                lp.scalarValue() += lps.scalarValue();
            }
        }

        //- create the surface modifier and move the surface points
        meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
        LongList<parMapperHelper> parallelBndNodes;

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(boundaryPoints, bpI)
        {
            labelledPointScalar& lps = preMapPositions[bpI];

            lps.coordinates() /= lps.scalarValue();

            const point& p = points[boundaryPoints[bpI]];

            //label patch, nearestTri;
            point pMap = p;
            scalar dSq;

            if( boundaryPointPatches[bpI].size() == 1 )
            {
                label nt;
                meshOctree_.findNearestSurfacePointInRegion
                (
                    pMap,
                    dSq,
                    nt,
                    boundaryPointPatches[bpI][0],
                    lps.coordinates()
                );
            }
            else
            {
                meshOctree_.findNearestPointToPatches
                (
                    pMap,
                    dSq,
                    lps.coordinates(),
                    boundaryPointPatches[bpI]
                );
            }

            const point newP = p + 0.5 * (pMap - p);

            surfaceModifier.moveBoundaryVertexNoUpdate(bpI, newP);

            if( bpAtProcsPtr && bpAtProcsPtr->sizeOfRow(bpI) )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                parallelBndNodes.append
                (
                    parMapperHelper
                    (
                        newP,
                        dSq,
                        bpI,
                        -1
                    )
                );
            }
        }

        //- make sure that the vertices at inter-processor boundaries
        //- are mapped onto the same location
        mapToSmallestDistance(parallelBndNodes);

        //- update the surface geometry of the
        surfaceModifier.updateGeometry();

        meshSurfaceOptimizer(surfaceEngine_, meshOctree_).untangleSurface();

        surfaceModifier.updateGeometry();
    }

    Info << "Finished smoothing mesh surface before mapping." << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
