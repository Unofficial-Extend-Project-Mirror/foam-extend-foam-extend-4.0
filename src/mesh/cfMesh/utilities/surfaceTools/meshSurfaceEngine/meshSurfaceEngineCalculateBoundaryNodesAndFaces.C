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
#include "boolList.H"
#include "helperFunctions.H"
#include "VRWGraphSMPModifier.H"
#include "labelledPoint.H"
#include "HashSet.H"

#include <map>
#include <set>

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceEngine::calculateBoundaryFaces() const
{
    if( mesh_.boundaries().size() != 0 )
    {
        const faceListPMG& faces = mesh_.faces();
        const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

        label nBoundaryFaces(0);
        if( activePatch_ < 0 )
        {
            //- take all patches
            forAll(boundaries, patchI)
                nBoundaryFaces += boundaries[patchI].patchSize();

            boundaryFacesPtr_ =
                new faceList::subList
                (
                    faces,
                    nBoundaryFaces,
                    boundaries[0].patchStart()
                );
        }
        else if( activePatch_ < boundaries.size() )
        {
            nBoundaryFaces = boundaries[activePatch_].patchSize();

            boundaryFacesPtr_ =
                new faceList::subList
                (
                    faces,
                    nBoundaryFaces,
                    boundaries[activePatch_].patchStart()
                );
        }
        else
        {
            FatalErrorIn
            (
                "void meshSurfaceEngine::calculateBoundaryFaces() const"
            ) << "Cannot select boundary faces. Invalid patch index "
              << activePatch_ << exit(FatalError);
        }

        reduce(nBoundaryFaces, sumOp<label>());
        Info << "Found " << nBoundaryFaces << " boundary faces " << endl;
    }
    else
    {
        FatalErrorIn
        (
            "void meshSurfaceEngine::calculateBoundaryFaces() const"
        ) << "Boundary faces are not at the end of the face list!"
            << exit(FatalError);
    }
}

void meshSurfaceEngine::calculateBoundaryOwners() const
{
    const labelList& owner = mesh_.owner();

    const faceList::subList& boundaryFaces = this->boundaryFaces();

    if( !boundaryFaceOwnersPtr_ )
    boundaryFaceOwnersPtr_ = new labelList(boundaryFaces.size());

    labelList& owners = *boundaryFaceOwnersPtr_;

    const label start = mesh_.boundaries()[0].patchStart();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(boundaryFaces, fI)
    owners[fI] = owner[start+fI];
}

void meshSurfaceEngine::calculateBoundaryNodes() const
{
    //- mark boundary points
    label pointI(0);
    if( !bppPtr_ )
        bppPtr_ = new labelList(mesh_.points().size(), -1);
    labelList& bp = *bppPtr_;

    const faceList::subList& boundaryFaces = this->boundaryFaces();

    boolList isBndPoint(bp.size(), false);

    # ifdef USE_OMP
    const label nThreads = 3 * omp_get_num_procs();
    # pragma omp parallel for num_threads(nThreads) schedule(static, 1)
    # endif
    forAll(boundaryFaces, bfI)
    {
        const face& bf = boundaryFaces[bfI];

        forAll(bf, pI)
            isBndPoint[bf[pI]] = true;
    }

    forAll(isBndPoint, pI)
    {
        if( isBndPoint[pI] )
            bp[pI] = pointI++;
    }

    if( Pstream::parRun() )
    {
        const faceListPMG& faces = mesh_.faces();
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();

        //- exchange information with processors
        //- this is need sometimes to find all nodes at the boundary
        bool found;
        do
        {
            found = false;

            //- send bnd nodes to other processor
            forAll(procBoundaries, patchI)
            {
                const label start = procBoundaries[patchI].patchStart();
                const label end = start + procBoundaries[patchI].patchSize();

                labelLongList dts;
                labelHashSet addedPoint;
                for(label faceI=start;faceI<end;++faceI)
                {
                    const face& f = faces[faceI];
                    forAll(f, pI)
                        if( (bp[f[pI]] != -1) && !addedPoint.found(f[pI]) )
                        {
                            addedPoint.insert(f[pI]);

                            //- data is sent as follows
                            //- 1. local face label in patch
                            //- 2. local node in face
                            dts.append(faceI-start);
                            dts.append((f.size() - pI) % f.size());
                        }
                }

                OPstream toOtherProc
                (
                    Pstream::commsTypes::blocking,
                    procBoundaries[patchI].neiProcNo(),
                    dts.byteSize()
                );
                toOtherProc << dts;
            }

            //- receive data and update positions if needed
            forAll(procBoundaries, patchI)
            {
                const label start = procBoundaries[patchI].patchStart();
                labelList receiveData;
                IPstream fromOtherProc
                (
                    Pstream::commsTypes::blocking,
                    procBoundaries[patchI].neiProcNo()
                );
                fromOtherProc >> receiveData;

                label counter(0);
                while( counter < receiveData.size() )
                {
                    const label fI = receiveData[counter++];
                    const label pI = receiveData[counter++];

                    if( bp[faces[start+fI][pI]] == -1 )
                    {
                        bp[faces[start+fI][pI]] = pointI++;
                        found = true;
                    }
                }
            }

            reduce(found, maxOp<bool>());

        } while( found );
    }

    if( !boundaryPointsPtr_ )
        boundaryPointsPtr_ = new labelList();

    labelList& boundaryPoints = *boundaryPointsPtr_;
    boundaryPoints.setSize(pointI);

    //- fill the boundaryPoints list
    # ifdef USE_OMP
    # pragma omp parallel for num_threads(nThreads) schedule(static, 1)
    # endif
    forAll(bp, bpI)
    {
        if( bp[bpI] != -1 )
            boundaryPoints[bp[bpI]] = bpI;
    }
}

void meshSurfaceEngine::calculateBoundaryFacePatches() const
{
    const faceList::subList& bFaces = this->boundaryFaces();
    boundaryFacePatchPtr_ = new labelList(bFaces.size());
    labelList& facePatch = *boundaryFacePatchPtr_;

    label faceI(0);
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    forAll(boundaries, patchI)
    {
        const label nFaces = boundaries[patchI].patchSize();
        for(label patchFaceI=0;patchFaceI<nFaces;++patchFaceI)
        {
            facePatch[faceI] = patchI;
            ++faceI;
        }
    }
}

void meshSurfaceEngine::calculatePointFaces() const
{
    //- fill pointFacesAddr
    if( !pointFacesPtr_ )
        pointFacesPtr_ = new VRWGraph();
    VRWGraph& pointFacesAddr = *pointFacesPtr_;

    if( !pointInFacePtr_ )
        pointInFacePtr_ = new VRWGraph();
    VRWGraph& pointInFaceAddr = *pointInFacePtr_;

    const labelList& bPoints = this->boundaryPoints();
    const faceList::subList& bFaces = this->boundaryFaces();

    //- create boundary points
    const labelList& bp = this->bp();

    labelLongList npf;

    # ifdef USE_OMP
    label nThreads = 3 * omp_get_num_procs();
    if( bPoints.size() < 1000 )
        nThreads = 1;
    # else
    const label nThreads(1);
    # endif

    label minRow(INT_MAX), maxRow(0);
    List<List<LongList<labelPair> > > dataForOtherThreads(nThreads);

    # ifdef USE_OMP
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI(0);
        # endif

        List<LongList<labelPair> >& dot = dataForOtherThreads[threadI];
        dot.setSize(nThreads);

        //- find min and max entry in the graph
        //- they are used for assigning ranges of values local for each process
        label localMinRow(minRow), localMaxRow(0);
        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];
            forAll(bf, pI)
            {
                const label bpI = bp[bf[pI]];
                localMaxRow = Foam::max(localMaxRow, bpI);
                localMinRow = Foam::min(localMinRow, bpI);
            }
        }

        ++localMaxRow;

        # ifdef USE_OMP
        # pragma omp critical
        # endif
        {
            minRow = Foam::min(minRow, localMinRow);
            minRow = Foam::max(minRow, 0);
            maxRow = Foam::max(maxRow, localMaxRow);

            npf.setSize(maxRow);
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        //- initialise appearances
        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        for(label i=0;i<maxRow;++i)
            npf[i] = 0;

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        const label range = (maxRow - minRow) / nThreads + 1;
        const label localMin = minRow + threadI * range;
        const label localMax = Foam::min(localMin + range, maxRow);

        //- find the number of appearances of each element in the original graph
        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            forAll(bf, pI)
            {
                const label bpI = bp[bf[pI]];

                const label threadNo = (bpI - minRow) / range;

                if( threadNo == threadI )
                {
                    ++npf[bpI];
                }
                else
                {
                    dot[threadNo].append(labelPair(bpI, bfI));
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        //- count the appearances which are not local to the processor
        for(label i=0;i<nThreads;++i)
        {
            const LongList<labelPair>& data =
                dataForOtherThreads[i][threadI];

            forAll(data, j)
                ++npf[data[j].first()];
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        //- allocate graph
        # ifdef USE_OMP
        # pragma omp master
        # endif
        {
            VRWGraphSMPModifier(pointFacesAddr).setSizeAndRowSize(npf);
            VRWGraphSMPModifier(pointInFaceAddr).setSizeAndRowSize(npf);
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        for(label i=localMin;i<localMax;++i)
            npf[i] = 0;

        //- start filling reverse addressing graph
        //- update data from processors with smaller labels
        for(label i=0;i<threadI;++i)
        {
            const LongList<labelPair>& data =
                dataForOtherThreads[i][threadI];

            forAll(data, j)
            {
                const label bpI = data[j].first();
                const label bfI = data[j].second();

                pointFacesAddr(bpI, npf[bpI]) = bfI;
                pointInFaceAddr(bpI, npf[bpI]) =
                    bFaces[bfI].which(bPoints[bpI]);

                ++npf[bpI];
            }
        }

        //- update data local to the processor
        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            forAll(bf, pI)
            {
                const label bpI = bp[bf[pI]];

                if( (bpI >= localMin) && (bpI < localMax) )
                {
                    pointInFaceAddr(bpI, npf[bpI]) = pI;
                    pointFacesAddr(bpI, npf[bpI]++) = bfI;
                }
            }
        }

        //- update data from the processors with higher labels
        for(label i=threadI+1;i<nThreads;++i)
        {
            const LongList<labelPair>& data =
                dataForOtherThreads[i][threadI];

            forAll(data, j)
            {
                const label bpI = data[j].first();
                const label bfI = data[j].second();

                pointFacesAddr(bpI, npf[bpI]) = bfI;
                pointInFaceAddr(bpI, npf[bpI]) =
                    bFaces[bfI].which(bPoints[bpI]);

                ++npf[bpI];
            }
        }
    }

    pointFacesAddr.setSize(bPoints.size());
    pointInFaceAddr.setSize(bPoints.size());
}

void meshSurfaceEngine::calculatePointPatches() const
{
    if( !pointPatchesPtr_ )
        pointPatchesPtr_ = new VRWGraph();
    VRWGraph& pPatches = *pointPatchesPtr_;

    const labelList& facePatch = boundaryFacePatches();
    const VRWGraph& pFaces = pointFaces();

    # ifdef USE_OMP
    const label nThreads = 3 * omp_get_num_procs();
    # endif

    labelList npPatches(pFaces.size());

    # ifdef USE_OMP
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(npPatches, i)
            npPatches[i] = 0;

        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(pFaces, bpI)
        {
            DynList<label> pf;
            forAllRow(pFaces, bpI, pfI)
                pf.appendIfNotIn(facePatch[pFaces(bpI, pfI)]);

            npPatches[bpI] = pf.size();
        }

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp master
        # endif
        VRWGraphSMPModifier(pPatches).setSizeAndRowSize(npPatches);

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp for schedule(static)
        # endif
        forAll(pFaces, bpI)
        {
            DynList<label> pf;
            forAllRow(pFaces, bpI, pfI)
                pf.appendIfNotIn(facePatch[pFaces(bpI, pfI)]);

            pPatches.setRow(bpI, pf);
        }
    }

    if( Pstream::parRun() )
    {
        const labelList& globalPointLabel = globalBoundaryPointLabel();
        const VRWGraph& bpAtProcs = this->bpAtProcs();
        const Map<label>& globalToLocal = globalToLocalBndPointAddressing();

        std::map<label, labelLongList> exchangeData;

        forAllConstIter(Map<label>, globalToLocal, iter)
        {
            const label bpI = iter();

            forAllRow(bpAtProcs, bpI, procI)
            {
                const label neiProc = bpAtProcs(bpI, procI);
                if( neiProc == Pstream::myProcNo() )
                    continue;

                if( exchangeData.find(neiProc) == exchangeData.end() )
                {
                    exchangeData.insert
                    (
                        std::make_pair(neiProc, labelLongList())
                    );
                }

                labelLongList& dataToSend = exchangeData[neiProc];

                //- prepare data which will be sent
                //- data is sent as follows
                //- 1. global point label
                //- 2. number of local patches for point
                //- 3. patch labels for a given point
                dataToSend.append(globalPointLabel[bpI]);
                dataToSend.append(pPatches.sizeOfRow(bpI));
                forAllRow(pPatches, bpI, patchI)
                    dataToSend.append(pPatches(bpI, patchI));
            }
        }

        //- exchange data with other processors
        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label bpI = globalToLocal[receivedData[counter++]];
            const label nPatches = receivedData[counter++];
            for(label i=0;i<nPatches;++i)
                pPatches.appendIfNotIn(bpI, receivedData[counter++]);
        }
    }
}

void meshSurfaceEngine::calculatePointPoints() const
{
    if( !pointPointsPtr_ )
        pointPointsPtr_ = new VRWGraph();

    VRWGraph& pointPoints = *pointPointsPtr_;

    const labelList& boundaryPoints = this->boundaryPoints();
    const faceList::subList& bFaces = this->boundaryFaces();
    const VRWGraph& pFaces = this->pointFaces();
    const labelList& bp = this->bp();

    # ifdef USE_OMP
    const label nThreads = 3 * omp_get_num_procs();
    # endif

    labelList npp(boundaryPoints.size());

    # ifdef USE_OMP
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(npp, i)
            npp[i] = 0;

        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(pFaces, bpI)
        {
            DynList<label> pPoints;

            forAllRow(pFaces, bpI, pfI)
            {
                const face& bf = bFaces[pFaces(bpI, pfI)];

                const label pos = bf.which(boundaryPoints[bpI]);

                pPoints.appendIfNotIn(bp[bf.nextLabel(pos)]);
                pPoints.appendIfNotIn(bp[bf.prevLabel(pos)]);
            }

            npp[bpI] = pPoints.size();
        }

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp master
        # endif
        VRWGraphSMPModifier(pointPoints).setSizeAndRowSize(npp);

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp for schedule(static)
        # endif
        forAll(pFaces, bpI)
        {
            DynList<label> pPoints;

            forAllRow(pFaces, bpI, pfI)
            {
                const face& bf = bFaces[pFaces(bpI, pfI)];

                const label pos = bf.which(boundaryPoints[bpI]);

                pPoints.appendIfNotIn(bp[bf.nextLabel(pos)]);
                pPoints.appendIfNotIn(bp[bf.prevLabel(pos)]);
            }

            pointPoints.setRow(bpI, pPoints);
        }
    }

    if( Pstream::parRun() )
    {
        //- this is needed to make the connection matrix symmetric
        //- on all processors. In some cases the points on a given processor
        //- may not be connected because of a single layer of faces on some
        //- other processor. P0, P0 | P1 | P0 P0
        const labelList& globalPointLabel = this->globalBoundaryPointLabel();
        const Map<label>& globalToLocal =
            this->globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = this->bpAtProcs();

        std::map<label, labelLongList> exchangeData;
        forAllConstIter(Map<label>, globalToLocal, iter)
        {
            const label bpI = iter();

            DynList<label> neiToSend;
            forAllRow(pointPoints, bpI, j)
            {
                const label bpJ = pointPoints(bpI, j);
                if( bpAtProcs.sizeOfRow(bpJ) != 0 )
                    neiToSend.append(bpJ);
            }

            forAllRow(bpAtProcs, bpI, procI)
            {
                const label neiProc = bpAtProcs(bpI, procI);
                if( neiProc == Pstream::myProcNo() )
                    continue;

                if( exchangeData.find(neiProc) == exchangeData.end() )
                    exchangeData.insert(std::make_pair(neiProc,labelLongList()));

                if( neiToSend.size() != 0 )
                {
                    labelLongList& dts = exchangeData[neiProc];
                    dts.append(globalPointLabel[bpI]);
                    dts.append(neiToSend.size());
                    forAll(neiToSend, i)
                        dts.append(globalPointLabel[neiToSend[i]]);
                }
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label bpI = globalToLocal[receivedData[counter++]];
            const label size = receivedData[counter++];
            for(label i=0;i<size;++i)
            {
                const label gpI = receivedData[counter++];
                if( globalToLocal.found(gpI) )
                    pointPoints.appendIfNotIn(bpI, globalToLocal[gpI]);
            }
        }
    }
}

void meshSurfaceEngine::calculatePointNormals() const
{
    const VRWGraph& pFaces = pointFaces();
    const vectorField& fNormals = faceNormals();

    pointNormalsPtr_ = new vectorField(pFaces.size());

    # ifdef USE_OMP
    # pragma omp parallel for if( pFaces.size() > 1000 ) schedule(dynamic, 50)
    # endif
    forAll(pFaces, pI)
    {
        vector normal(vector::zero);

        forAllRow(pFaces, pI, pfI)
            normal += fNormals[pFaces(pI, pfI)];

        const scalar d = mag(normal);
        if( d > VSMALL )
        {
            normal /= d;
        }
        else
        {
            normal = vector::zero;
        }

        (*pointNormalsPtr_)[pI] = normal;
    }

    updatePointNormalsAtProcBoundaries();
}

void meshSurfaceEngine::calculateFaceNormals() const
{
    const faceList::subList& bFaces = this->boundaryFaces();
    const pointFieldPMG& points = mesh_.points();

    faceNormalsPtr_ = new vectorField(bFaces.size());

    # ifdef USE_OMP
    # pragma omp parallel for if( bFaces.size() > 1000 )
    # endif
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        faceNormalsPtr_->operator[](bfI) = bf.normal(points);
    }
}

void meshSurfaceEngine::calculateFaceCentres() const
{
    const faceList::subList& bFaces = this->boundaryFaces();
    const pointFieldPMG& points = mesh_.points();

    faceCentresPtr_ = new vectorField(bFaces.size());

    # ifdef USE_OMP
    # pragma omp parallel for if( bFaces.size() > 1000 )
    # endif
    forAll(bFaces, bfI)
        faceCentresPtr_->operator[](bfI) = bFaces[bfI].centre(points);
}

void meshSurfaceEngine::updatePointNormalsAtProcBoundaries() const
{
    if( !Pstream::parRun() )
        return;

    const VRWGraph& pFaces = pointFaces();
    const vectorField& fNormals = faceNormals();
    const labelList& globalPointLabel = this->globalBoundaryPointLabel();
    const Map<label>& globalToLocal =
        this->globalToLocalBndPointAddressing();
    const VRWGraph& bpAtProcs = this->bpAtProcs();

    vectorField& pNormals = *pointNormalsPtr_;

    //- create data which will be sent to other processors
    std::map<label, LongList<labelledPoint> > exchangeData;

    forAllConstIter(Map<label>, globalToLocal, iter)
    {
        const label bpI = iter();

        vector& n = pNormals[bpI];
        n = vector::zero;

        forAllRow(pFaces, bpI, pfI)
            n += fNormals[pFaces(bpI, pfI)];

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;
            if( exchangeData.find(neiProc) == exchangeData.end() )
            {
                exchangeData.insert
                (
                    std::make_pair(neiProc, LongList<labelledPoint>())
                );
            }

            exchangeData[neiProc].append
            (
                labelledPoint(globalPointLabel[bpI], n)
            );
        }
    }

    //- exchange data with other procs
    LongList<labelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    forAll(receivedData, i)
    {
        const label bpI = globalToLocal[receivedData[i].pointLabel()];
        pNormals[bpI] += receivedData[i].coordinates();
    }

    //- normalize vectors
    # ifdef USE_OMP
    # pragma omp parallel for if( bpAtProcs.size() > 1000 ) \
    schedule(guided)
    # endif
    forAll(bpAtProcs, bpI)
    {
        if( bpAtProcs.sizeOfRow(bpI) == 0 )
            continue;

        vector normal = pNormals[bpI];
        const scalar d = mag(normal);
        if( d > VSMALL )
        {
            normal /= d;
        }
        else
        {
            normal = vector::zero;
        }

        pNormals[bpI] = normal;
    }
}

void meshSurfaceEngine::calculateEdgesAndAddressing() const
{
    const VRWGraph& pFaces = pointFaces();
    const faceList::subList& bFaces = boundaryFaces();
    const labelList& bPoints = boundaryPoints();
    const labelList& bp = this->bp();

    edgesPtr_ = new edgeList();
    edgeList& edges = *edgesPtr_;

    bpEdgesPtr_ = new VRWGraph();
    VRWGraph& bpEdges = *bpEdgesPtr_;

    # ifdef USE_OMP
    label nThreads = 3 * omp_get_num_procs();
    if( pFaces.size() < 1000 )
        nThreads = 1;
    # else
    const label nThreads(1);
    # endif

    labelList nEdgesForThread(nThreads);

    label edgeI(0);

    # ifdef USE_OMP
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        edgeLongList edgesHelper;

        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(pFaces, bpI)
        {
            std::set<std::pair<label, label> > edgesAtPoint;

            forAllRow(pFaces, bpI, pfI)
            {
                const label bfI = pFaces(bpI, pfI);
                const face& bf = bFaces[bfI];

                const label pos = bf.which(bPoints[bpI]);

                if( bp[bf.nextLabel(pos)] >= bpI )
                {
                    edgesAtPoint.insert
                    (
                        std::make_pair(bf[pos], bf.nextLabel(pos))
                    );
                }
                if( bp[bf.prevLabel(pos)] >= bpI )
                {
                    edgesAtPoint.insert
                    (
                        std::make_pair(bf[pos], bf.prevLabel(pos))
                    );
                }
            }

            std::set<std::pair<label, label> >::const_iterator it;
            for(it=edgesAtPoint.begin();it!=edgesAtPoint.end();++it)
                edgesHelper.append(edge(it->first, it->second));
        }

        //- this enables other threads to see the number of edges
        //- generated by each thread
        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI(0);
        # endif
        nEdgesForThread[threadI] = edgesHelper.size();

        # ifdef USE_OMP
        # pragma omp critical
        # endif
        edgeI += edgesHelper.size();

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp master
        # endif
        edgesPtr_->setSize(edgeI);

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        //- find the starting position of the edges generated by this thread
        //- in the global list of edges
        label localStart(0);
        for(label i=0;i<threadI;++i)
            localStart += nEdgesForThread[i];

        //- store edges into the global list
        forAll(edgesHelper, i)
            edgesPtr_->operator[](localStart+i) = edgesHelper[i];
    }

    //- set the bpEdges
    VRWGraphSMPModifier(bpEdges).reverseAddressing(bp, edges);
    bpEdges.setSize(pFaces.size());

    if( !Pstream::parRun() )
        return;

    bool addEdges;
    do
    {
        addEdges = false;

        //- mark boundary edges for processors which do not contain
        //- boundary faces. This procedure is needed to identify boundary
        //- edges which are not part of any boundary face on their processor
        const faceListPMG& faces = mesh_.faces();
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();

        //- send boundary edges to neighbour processors
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label end = start + procBoundaries[patchI].patchSize();

            labelLongList dts;
            for(label faceI=start;faceI<end;++faceI)
            {
                const face& f = faces[faceI];
                forAll(f, eI)
                {
                    const edge e = f.faceEdge(eI);
                    const label s = bp[e.start()];
                    if( s < 0 )
                        continue;

                    forAllRow(bpEdges, s, peI)
                        if( edges[bpEdges(s, peI)] == e )
                        {
                            dts.append(faceI-start);
                            dts.append((f.size()-1-eI)%f.size());
                            break;
                        }
                }
            }

            OPstream toOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                dts.byteSize()
            );
            toOtherProc << dts;
        }

        //- receive data from other processors. Mark edges which are not yet
        //- marked as boundary edges
        forAll(procBoundaries, patchI)
        {
            labelList receivedEdges;
            IPstream fromOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo()
            );
            fromOtherProc >> receivedEdges;

            const label start = procBoundaries[patchI].patchStart();
            label nReceivedEdges(0);
            while( nReceivedEdges < receivedEdges.size() )
            {
                const face& f = faces[start+receivedEdges[nReceivedEdges++]];
                const label eI = receivedEdges[nReceivedEdges++];

                const edge e = f.faceEdge(eI);
                const label s = bp[e.start()];

                bool found(false);
                forAllRow(bpEdges, s, peI)
                    if( edges[bpEdges(s, peI)] == e )
                    {
                        found = true;
                        break;
                    }

                if( !found )
                {
                    //- create a new edge
                    addEdges = true;
                    edges.newElmt(edgeI) = e;

                    bpEdges.append(bp[e.start()], edgeI);
                    bpEdges.append(bp[e.end()], edgeI);
                    ++edgeI;
                }
            }
        }

        reduce(addEdges, maxOp<bool>());
    } while( addEdges );

    edges.setSize(edgeI);
}

void meshSurfaceEngine::calculateFaceEdgesAddressing() const
{
    const faceList::subList& bFaces = this->boundaryFaces();
    const labelList& bp = this->bp();
    const edgeList& edges = this->edges();
    const VRWGraph& bpEdges = this->boundaryPointEdges();

    faceEdgesPtr_ = new VRWGraph(bFaces.size());
    VRWGraph& faceEdges = *faceEdgesPtr_;

    labelList nfe(bFaces.size());

    # ifdef USE_OMP
    const label nThreads = 3 * omp_get_num_procs();

    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(bFaces, bfI)
            nfe[bfI] = bFaces[bfI].size();

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp master
        {
        # endif

        VRWGraphSMPModifier(faceEdges).setSizeAndRowSize(nfe);

        # ifdef USE_OMP
        }

        # pragma omp barrier

        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(faceEdges, bfI)
        {
            const face& bf = bFaces[bfI];

            forAll(bf, eI)
            {
                const edge e = bf.faceEdge(eI);

                const label bps = bp[e.start()];

                forAllRow(bpEdges, bps, peI)
                {
                    const label beI = bpEdges(bps, peI);
                    const edge& ee = edges[beI];

                    if( e == ee )
                    {
                        faceEdges(bfI, eI) = beI;
                        break;
                    }
                }
            }
        }
    }
}

void meshSurfaceEngine::calculateEdgeFacesAddressing() const
{
    const faceList::subList& bFaces = this->boundaryFaces();
    const VRWGraph& pointFaces = this->pointFaces();
    const edgeList& edges = this->edges();
    const labelList& bp = this->bp();

    edgeFacesPtr_ = new VRWGraph();
    VRWGraph& edgeFaces = *edgeFacesPtr_;

    labelList nef(edges.size());

    # ifdef USE_OMP
    const label nThreads = 3 * omp_get_num_procs();

    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(nef, edgeI)
            nef[edgeI] = 0;

        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(edges, edgeI)
        {
            const edge& ee = edges[edgeI];
            const label bpI = bp[ee.start()];

            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);

                const face& bf = bFaces[bfI];

                forAll(bf, eI)
                {
                    if( bf.faceEdge(eI) == ee )
                    {
                        ++nef[edgeI];
                        break;
                    }
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp master
        # endif
        VRWGraphSMPModifier(edgeFaces).setSizeAndRowSize(nef);

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp for schedule(static)
        # endif
        forAll(edges, edgeI)
        {
            const edge& ee = edges[edgeI];
            const label bpI = bp[ee.start()];

            //- find boundary faces attached to this edge
            DynList<label> eFaces;
            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);

                const face& bf = bFaces[bfI];

                forAll(bf, eI)
                {
                    if( bf.faceEdge(eI) == ee )
                    {
                        eFaces.append(bfI);
                        break;
                    }
                }
            }

            //- the face that owns the edge shall be the first one in the list
            // TODO: find out whether this will be necessary
            if( eFaces.size() == 2 )
            {
                const face& bf = bFaces[eFaces[1]];

                const label pos = bf.which(ee.start());

                if( bf.nextLabel(pos) == ee.end() )
                {
                    //- this face shall be the first one in the list
                    const label helper = eFaces[0];
                    eFaces[0] = eFaces[1];
                    eFaces[1] = helper;
                }
            }

            edgeFaces.setRow(edgeI, eFaces);
        }
    }
}

void meshSurfaceEngine::calculateEdgePatchesAddressing() const
{
    edgePatchesPtr_ = new VRWGraph();
    VRWGraph& edgePatches = *edgePatchesPtr_;

    const VRWGraph& edgeFaces = this->edgeFaces();
    const labelList& facePatch = this->boundaryFacePatches();

    edgePatches.setSize(edgeFaces.size());

    forAll(edgeFaces, eI)
    {
        DynList<label> ePatches;

        forAllRow(edgeFaces, eI, i)
        {
            const label patchI = facePatch[edgeFaces(eI, i)];

            ePatches.appendIfNotIn(patchI);
        }

        edgePatches.setRow(eI, ePatches);
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal = globalToLocalBndEdgeAddressing();
        const Map<label>& otherPatch = this->otherEdgeFacePatch();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label beI = it();

            edgePatches.appendIfNotIn(beI, otherPatch[beI]);
        }
    }
}

void meshSurfaceEngine::calculateFaceFacesAddressing() const
{
    const VRWGraph& edgeFaces = this->edgeFaces();

    const faceList::subList& bFaces = boundaryFaces();
    faceFacesPtr_ = new VRWGraph(bFaces.size());
    VRWGraph& faceFaces = *faceFacesPtr_;

    forAll(bFaces, bfI)
        faceFaces.setRowSize(bfI, bFaces[bfI].size());

    labelList nAppearances(bFaces.size(), 0);

    forAll(edgeFaces, efI)
    {
        if( edgeFaces.sizeOfRow(efI) == 2 )
        {
            const label f0 = edgeFaces(efI, 0);
            const label f1 = edgeFaces(efI, 1);

            faceFaces(f0, nAppearances[f0]++) = f1;
            faceFaces(f1, nAppearances[f1]++) = f0;
        }
        else if( Pstream::parRun() && (edgeFaces.sizeOfRow(efI) == 1) )
        {
            const label f0 = edgeFaces(efI, 0);
            faceFaces(f0, nAppearances[f0]++) = -1;
        }
        else if( Pstream::parRun() && (edgeFaces.sizeOfRow(efI) != 0 ) )
        {
            FatalErrorIn
            (
                "void meshSurfaceEngine::calculateFaceFacesAddressing() const"
            ) << "The surface of the mesh is invalid!"
                << " The number of faces containing edge " << efI
                << " is " << edgeFaces.sizeOfRow(efI)
                << " Cannot continue" << exit(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
