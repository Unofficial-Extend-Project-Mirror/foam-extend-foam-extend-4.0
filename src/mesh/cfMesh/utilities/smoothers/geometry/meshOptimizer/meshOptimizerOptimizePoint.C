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
#include "meshOptimizer.H"
#include "polyMeshGenAddressing.H"
#include "meshSurfaceEngine.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOptimizer::laplaceSmoother::laplacian
(
    const labelLongList& smoothPoints,
    const label nIterations
)
{
    const VRWGraph& pPoints = mesh_.addressingData().pointPoints();
    pointFieldPMG& points = mesh_.points();

    for(label iterationI=0;iterationI<nIterations;++iterationI)
    {
        labelLongList procPoints;

        forAll(smoothPoints, i)
        {
            const label pointI = smoothPoints[i];

            if( vertexLocation_[pointI] & LOCKED )
                continue;

            if( vertexLocation_[pointI] & PARALLELBOUNDARY )
            {
                procPoints.append(pointI);

                continue;
            }

            vector newP(vector::zero);

            const label nPointPoints = pPoints.sizeOfRow(pointI);

            if( nPointPoints == 0 )
                return;

            for(label pI=0;pI<nPointPoints;++pI)
                newP += points[pPoints(pointI, pI)];

            newP /= pPoints.sizeOfRow(pointI);
            points[pointI] = newP;
        }

        laplacianParallel(procPoints, false);
    }

    updateMeshGeometry(smoothPoints);
}

void meshOptimizer::laplaceSmoother::laplacianSurface
(
    const labelLongList& smoothPoints,
    const label nIterations
)
{
    const VRWGraph& pPoints = mesh_.addressingData().pointPoints();
    pointFieldPMG& points = mesh_.points();

    for(label iterationI=0;iterationI<nIterations;++iterationI)
    {
        labelLongList procPoints;

        forAll(smoothPoints, i)
        {
            const label pointI = smoothPoints[i];

            if( vertexLocation_[pointI] & LOCKED )
                continue;

            if( vertexLocation_[pointI] & PARALLELBOUNDARY )
            {
                procPoints.append(pointI);

                continue;
            }

            vector newP(vector::zero);

            label counter(0);
            forAllRow(pPoints, pointI, pI)
            {
                const label pLabel = pPoints(pointI, pI);
                if( vertexLocation_[pLabel] & INSIDE )
                    continue;

                newP += points[pLabel];
                ++counter;
            }

            if( counter != 0 )
            {
                newP /= counter;
                points[pointI] = newP;
            }
        }

        laplacianParallel(smoothPoints, true);
    }

    updateMeshGeometry(smoothPoints);
}

void meshOptimizer::laplaceSmoother::laplacianPC
(
    const labelLongList& smoothPoints,
    const label nIterations
)
{
    const VRWGraph& pointCells = mesh_.addressingData().pointCells();
    const vectorField& centres = mesh_.addressingData().cellCentres();
    pointFieldPMG& points = mesh_.points();

    for(label iterationI=0;iterationI<nIterations;++iterationI)
    {
        labelLongList procPoints;

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 20)
        # endif
        forAll(smoothPoints, i)
        {
            const label pointI = smoothPoints[i];

            if( vertexLocation_[pointI] & LOCKED )
                continue;

            if( pointCells.sizeOfRow(pointI) == 0 )
                continue;

            if( vertexLocation_[pointI] & PARALLELBOUNDARY )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                procPoints.append(pointI);

                continue;
            }

            point newP(vector::zero);
            forAllRow(pointCells, pointI, pcI)
                newP += centres[pointCells(pointI, pcI)];

            newP /= pointCells.sizeOfRow(pointI);

            points[pointI] = newP;
        }

        laplacianPCParallel(procPoints);

        updateMeshGeometry(smoothPoints);
    }
}

void meshOptimizer::laplaceSmoother::laplacianWPC
(
    const labelLongList& smoothPoints,
    const label nIterations
)
{
    const VRWGraph& pointCells = mesh_.addressingData().pointCells();
    const vectorField& centres = mesh_.addressingData().cellCentres();
    const scalarField& volumes = mesh_.addressingData().cellVolumes();

    pointFieldPMG& points = mesh_.points();

    for(label iterationI=0;iterationI<nIterations;++iterationI)
    {
        labelLongList procPoints;

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 20)
        # endif
        forAll(smoothPoints, i)
        {
            const label pointI = smoothPoints[i];

            if( vertexLocation_[pointI] & LOCKED )
                continue;

            if( pointCells.sizeOfRow(pointI) == 0 )
                continue;

            if( vertexLocation_[pointI] & PARALLELBOUNDARY )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                procPoints.append(pointI);

                continue;
            }

            point newP(vector::zero);
            scalar sumWeights(0.0);
            forAllRow(pointCells, pointI, pcI)
            {
                const label cellI = pointCells(pointI, pcI);
                const scalar w = Foam::max(volumes[cellI], VSMALL);
                newP += w * centres[cellI];
                sumWeights += w;
            }

            newP /= sumWeights;
            points[pointI] = newP;
        }

        laplacianWPCParallel(procPoints);

        updateMeshGeometry(smoothPoints);
    }
}

void meshOptimizer::laplaceSmoother::updateMeshGeometry
(
    const labelLongList& smoothPoints
)
{
    const cellListPMG& cells = mesh_.cells();
    const VRWGraph& pointCells = mesh_.addressingData().pointCells();

    boolList chF(mesh_.faces().size(), false);

    # ifdef USE_OMP
    # pragma omp parallel for if( smoothPoints.size() > 100 ) \
    schedule(dynamic, 20)
    # endif
    forAll(smoothPoints, i)
    {
        const label pointI = smoothPoints[i];

        if( vertexLocation_[pointI] & LOCKED )
            continue;

        forAllRow(pointCells, pointI, pcI)
        {
            const cell& c = cells[pointCells(pointI, pcI)];

            forAll(c, fI)
                chF[c[fI]] = true;
        }
    }

    //- make sure that neighbouring processors get the same information
    const PtrList<processorBoundaryPatch>& pBnd = mesh_.procBoundaries();
    forAll(pBnd, patchI)
    {
        const label start = pBnd[patchI].patchStart();
        const label size = pBnd[patchI].patchSize();

        labelLongList sendData;
        for(label faceI=0;faceI<size;++faceI)
        {
            if( chF[start+faceI] )
                sendData.append(faceI);
        }

        OPstream toOtherProc
        (
            Pstream::blocking,
            pBnd[patchI].neiProcNo(),
            sendData.byteSize()
        );

        toOtherProc << sendData;
    }

    forAll(pBnd, patchI)
    {
        labelList receivedData;

        IPstream fromOtherProc
        (
            Pstream::blocking,
            pBnd[patchI].neiProcNo()
        );

        fromOtherProc >> receivedData;

        const label start = pBnd[patchI].patchStart();
        forAll(receivedData, i)
            chF[start+receivedData[i]] = true;
    }

    //- update geometry information
    const_cast<polyMeshGenAddressing&>
    (
        mesh_.addressingData()
    ).updateGeometry(chF);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Constructor of laplaceSmoother

meshOptimizer::laplaceSmoother::laplaceSmoother
(
    polyMeshGen& mesh,
    const List<direction>& vertexLocation
)
:
    mesh_(mesh),
    vertexLocation_(vertexLocation)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

meshOptimizer::laplaceSmoother::~laplaceSmoother()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Member Functions

void meshOptimizer::laplaceSmoother::optimizeLaplacian(const label nIterations)
{
    labelLongList smoothPoints;

    forAll(vertexLocation_, pointI)
    {
        if( vertexLocation_[pointI] & INSIDE )
            smoothPoints.append(pointI);
    }

    laplacian(smoothPoints, nIterations);
}

void meshOptimizer::laplaceSmoother::optimizeLaplacian
(
    const labelHashSet& badFaces,
    const label nIterations
)
{
    FatalError << "Not implemented " << exit(FatalError);
}

void meshOptimizer::laplaceSmoother::optimizeSurfaceLaplacian
(
    const labelHashSet& badFaces,
    const label nIterations
)
{
    FatalError << "Not implemented " << exit(FatalError);
}

void meshOptimizer::laplaceSmoother::optimizeLaplacianPC
(
    const label nIterations
)
{
    labelLongList smoothPoints;

    forAll(vertexLocation_, pointI)
    {
        if( vertexLocation_[pointI] & INSIDE )
            smoothPoints.append(pointI);
    }

    laplacianPC(smoothPoints, nIterations);
}

void meshOptimizer::laplaceSmoother::optimizeLaplacianPC
(
    const labelHashSet& badFaces,
    const label nIterations
)
{
    FatalError << "Not implemented " << exit(FatalError);
}

void meshOptimizer::laplaceSmoother::optimizeLaplacianWPC
(
    const label nIterations
)
{
    labelLongList smoothPoints;

    forAll(vertexLocation_, pointI)
    {
        if( vertexLocation_[pointI] & INSIDE )
            smoothPoints.append(pointI);
    }

    laplacianWPC(smoothPoints, nIterations);
}

void meshOptimizer::laplaceSmoother::optimizeLaplacianWPC
(
    const labelHashSet& badFaces,
    const label nIterations
)
{
    FatalError << "Not implemented " << exit(FatalError);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
