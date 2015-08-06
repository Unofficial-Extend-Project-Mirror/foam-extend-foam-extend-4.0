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

#include "error.H"
#include "polyMeshGenModifier.H"
#include "edgeExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfaceOptimizer.H"
#include "meshOctree.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "helperFunctions.H"
#include "DynList.H"
#include "labelPair.H"
#include "labelledScalar.H"
#include "labelledPoint.H"
#include "refLabelledPoint.H"
#include "HashSet.H"
#include "triSurfacePartitioner.H"
#include "triSurfaceClassifyEdges.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGEdgeExtractor

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

void edgeExtractor::faceEvaluator::calculateNeiPatchesParallel()
{
    otherFacePatch_.clear();

    const labelList& fPatches = extractor_.facePatch_;

    if( Pstream::parRun() )
    {
        const meshSurfaceEngine& mse = extractor_.surfaceEngine();
        const VRWGraph& edgeFaces = mse.edgeFaces();
        const Map<label>& otherProc = mse.otherEdgeFaceAtProc();
        const Map<label>& globalToLocal =
            mse.globalToLocalBndEdgeAddressing();

        //- create communication matrix
        std::map<label, labelLongList> exchangeData;
        const DynList<label>& neiProcs = mse.beNeiProcs();
        forAll(neiProcs, procI)
            exchangeData.insert
            (
                std::make_pair(neiProcs[procI], labelLongList())
            );

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label beI = it();

            if( edgeFaces.sizeOfRow(beI) == 1 )
            {
                labelLongList& dts = exchangeData[otherProc[beI]];
                //- send data as follows:
                //- 1. global edge label
                //- 2. patch of the attached boundary face
                dts.append(it.key());
                dts.append(fPatches[edgeFaces(beI, 0)]);
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label beI = globalToLocal[receivedData[counter++]];
            const label fPatch = receivedData[counter++];

            otherFacePatch_.insert(beI, fPatch);
        }
    }
}

void edgeExtractor::faceEvaluator::calculateNeiPatchesParallelNewPatches()
{
    if( newOtherFacePatchPtr_ )
        return;

    if( !newBoundaryPatchesPtr_ )
        FatalErrorIn
        (
            "void edgeExtractor::faceEvaluator::"
            "calculateNeiPatchesParallelNewPatches()"
        ) << "newBoundaryPatchesPtr_ are NULL" << exit(FatalError);

    newOtherFacePatchPtr_ = new Map<label>();
    Map<label>& otherFacePatch = *newOtherFacePatchPtr_;

    const labelList& fPatches = *newBoundaryPatchesPtr_;

    if( Pstream::parRun() )
    {
        const meshSurfaceEngine& mse = extractor_.surfaceEngine();
        const VRWGraph& edgeFaces = mse.edgeFaces();
        const Map<label>& otherProc = mse.otherEdgeFaceAtProc();
        const Map<label>& globalToLocal =
            mse.globalToLocalBndEdgeAddressing();

        //- create communication matrix
        std::map<label, labelLongList> exchangeData;
        const DynList<label>& neiProcs = mse.beNeiProcs();
        forAll(neiProcs, procI)
            exchangeData.insert
            (
                std::make_pair(neiProcs[procI], labelLongList())
            );

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label beI = it();

            if( edgeFaces.sizeOfRow(beI) == 1 )
            {
                labelLongList& dts = exchangeData[otherProc[beI]];
                //- send data as follows:
                //- 1. global edge label
                //- 2. patch of the attached boundary face
                dts.append(it.key());
                dts.append(fPatches[edgeFaces(beI, 0)]);
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label beI = globalToLocal[receivedData[counter++]];
            const label fPatch = receivedData[counter++];

            otherFacePatch.insert(beI, fPatch);
        }
    }
}

void edgeExtractor::faceEvaluator::neiFacesOverEdges
(
    const label bfI,
    DynList<label>& neiFaces
) const
{
    const meshSurfaceEngine& mse = extractor_.surfaceEngine();

    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    neiFaces.setSize(faceEdges.sizeOfRow(bfI));

    forAllRow(faceEdges, bfI, feI)
    {
        const label beI = faceEdges(bfI, feI);

        if( edgeFaces.sizeOfRow(beI) == 2 )
        {
            neiFaces[feI] = edgeFaces(beI, 0);
            if( neiFaces[feI] == bfI )
                neiFaces[feI] = edgeFaces(beI, 1);
        }
        else
        {
            neiFaces[feI] = -1;
        }
    }
}

void edgeExtractor::faceEvaluator::neiFacesProcs
(
    const label bfI,
    DynList<label>& neiProcs
) const
{
    const meshSurfaceEngine& mse = extractor_.surfaceEngine();

    const VRWGraph& faceEdges = mse.faceEdges();

    neiProcs.setSize(faceEdges.sizeOfRow(bfI));
    neiProcs = Pstream::myProcNo();

    if( Pstream::parRun() )
    {
        const Map<label>& otherFaceProc = mse.otherEdgeFaceAtProc();

        forAllRow(faceEdges, bfI, feI)
        {
            const label beI = faceEdges(bfI, feI);

            const Map<label>::const_iterator it = otherFaceProc.find(beI);

            if( it != otherFaceProc.end() )
                neiProcs[feI] = it();
        }
    }
}

void edgeExtractor::faceEvaluator::neiPatchesOverEdges
(
    const label bfI,
    const labelList& fPatches,
    const Map<label>& otherFacePatch,
    DynList<label> &neiPatches
) const
{
    const meshSurfaceEngine& mse = extractor_.surfaceEngine();

    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    neiPatches.setSize(faceEdges.sizeOfRow(bfI));

    forAllRow(faceEdges, bfI, feI)
    {
        const label beI = faceEdges(bfI, feI);

        if( edgeFaces.sizeOfRow(beI) == 2 )
        {
            label nei = edgeFaces(beI, 0);
            if( nei == bfI )
                nei = edgeFaces(beI, 1);

            neiPatches[feI] = fPatches[nei];
        }
        else if( Pstream::parRun() && (edgeFaces.sizeOfRow(beI) == 1) )
        {
            neiPatches[feI] = otherFacePatch[beI];
        }
    }
}

label edgeExtractor::faceEvaluator::bestPatchTopological
(
    const DynList<label>& neiPatches,
    const label currentPatch
)
{
    //- find indices of all neighbour patches
    DynList<label> allNeiPatches;
    forAll(neiPatches, i)
        allNeiPatches.appendIfNotIn(neiPatches[i]);

    if( (allNeiPatches.size() == 1) && (allNeiPatches[0] == currentPatch) )
        return currentPatch;

    //- counter the number of neighbours in a patch
    Map<label> nNeiInPatch(allNeiPatches.size());
    forAll(allNeiPatches, i)
        nNeiInPatch.insert(allNeiPatches[i], 0);
    forAll(neiPatches, eI)
        ++nNeiInPatch[neiPatches[eI]];

    label newPatch = -1;
    label nNeiEdges(0);
    forAllConstIter(Map<label>, nNeiInPatch, it)
    {
        if( it() > nNeiEdges )
        {
            newPatch = it.key();
            nNeiEdges = it();
        }
        else if
        (
            (it() == nNeiEdges) && (it.key() == currentPatch)
        )
        {
            newPatch = it.key();
        }
    }

    //- do not swap if the situation allows for more than one edge
    //- shared with faces in other patches than the dominant one
    if( nNeiEdges < (neiPatches.size() - 1) )
        return currentPatch;

    //- do not swap if the best face is in the current patch
    if( (newPatch < 0) || (newPatch == currentPatch) )
        return currentPatch;

    //- check whether the edges shared ith the neighbour patch form
    //- a singly linked chain
    DynList<bool> sharedEdge;
    sharedEdge.setSize(neiPatches.size());
    sharedEdge = false;

    forAll(neiPatches, eI)
        if( neiPatches[eI] == newPatch )
            sharedEdge[eI] = true;

    if( help::areElementsInChain(sharedEdge) )
        return newPatch;

    return currentPatch;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

edgeExtractor::faceEvaluator::faceEvaluator(const edgeExtractor& ee)
:
    extractor_(ee),
    otherFacePatch_(),
    newBoundaryPatchesPtr_(NULL),
    newOtherFacePatchPtr_(NULL)
{
    if( Pstream::parRun() )
        calculateNeiPatchesParallel();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

edgeExtractor::faceEvaluator::~faceEvaluator()
{
    deleteDemandDrivenData(newOtherFacePatchPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

void edgeExtractor::faceEvaluator::setNewBoundaryPatches
(
    const labelList& newBoudaryPatches
)
{
    newBoundaryPatchesPtr_ = &newBoudaryPatches;

    if( Pstream::parRun() )
        calculateNeiPatchesParallelNewPatches();
}

void edgeExtractor::faceEvaluator::neiPatchesOverEdges
(
    const label bfI,
    DynList<label>& neiPatches
) const
{
    neiPatchesOverEdges
    (
        bfI,
        extractor_.facePatch_,
        otherFacePatch_,
        neiPatches
    );
}

label edgeExtractor::faceEvaluator::bestPatchTopological(const label bfI) const
{
    //- get neighbour patches over edges
    DynList<label> neiPatches;
    neiPatchesOverEdges
    (
        bfI,
        extractor_.facePatch_,
        otherFacePatch_,
        neiPatches
    );

    return bestPatchTopological(neiPatches, extractor_.facePatch_[bfI]);
}

label edgeExtractor::faceEvaluator::bestPatchAfterModification
(
    const label bfI
) const
{
    const label patchI = newBoundaryPatchesPtr_->operator[](bfI);

    if( patchI != extractor_.facePatch_[bfI] )
    {
        DynList<label> newNeiPatches, oldNeiPatches;
        neiPatchesOverEdges
        (
            bfI,
            *newBoundaryPatchesPtr_,
            *newOtherFacePatchPtr_,
            newNeiPatches
        );

        neiPatchesOverEdges
        (
            bfI,
            extractor_.facePatch_,
            otherFacePatch_,
            oldNeiPatches
        );

        DynList<label> neiFaces, neiProcs;
        neiFacesOverEdges(bfI, neiFaces);
        neiFacesProcs(bfI, neiProcs);

        forAll(neiFaces, eI)
        {
            const label origPatchI = oldNeiPatches[eI];
            const label newPatchI = newNeiPatches[eI];

            if( neiFaces[eI] > bfI )
            {
                if( newPatchI != origPatchI )
                    newNeiPatches[eI] = origPatchI;
            }
            else if( neiFaces[eI] < 0 )
            {
                if( neiProcs[eI] > Pstream::myProcNo() )
                    newNeiPatches[eI] = origPatchI;
            }
        }

        return bestPatchTopological(newNeiPatches, patchI);
    }

    return patchI;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

void edgeExtractor::cornerEvaluator::createParallelAddressing()
{
    const labelHashSet& corners = partitioner_.corners();

    const labelList& facePatch = extractor_.facePatch_;

    typedef Map<label> mapType;

    const meshSurfaceEngine& mse = extractor_.surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const pointFieldPMG& points = mse.points();
    const labelList& bp = mse.bp();
    const VRWGraph& pointFaces = mse.pointFaces();
    const labelList& globalLabel = mse.globalBoundaryPointLabel();
    const mapType& globalToLocal = mse.globalToLocalBndPointAddressing();
    const VRWGraph& bpAtProcs = mse.bpAtProcs();
    const DynList<label>& neiProcs = mse.bpNeiProcs();

    typedef std::map<label, LongList<labelledPoint> > exchangeMapType;
    exchangeMapType exchangeData;
    forAll(neiProcs, i)
        exchangeData[neiProcs[i]].clear();

    faceMap_.clear();
    faceAtProc_.clear();
    facePatches_.clear();

    forAllConstIter(mapType, globalToLocal, it)
    {
        const label bpI = it();

        if( corners.found(bpI) )
        {
            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                LongList<labelledPoint>& dts = exchangeData[neiProc];

                //- data is send in the ollowing form
                //- 1. globsl label of the current point
                //- 2. number of faces attached to the point
                //- 3. for each face send its patch, number of points and point
                dts.append(labelledPoint(it.key(), vector::zero));
                dts.append
                (
                    labelledPoint(pointFaces.sizeOfRow(bpI), vector::zero)
                );
                forAllRow(pointFaces, bpI, i)
                {
                    const face& bf = bFaces[pointFaces(bpI, i)];
                    dts.append
                    (
                        labelledPoint
                        (
                            facePatch[pointFaces(bpI, i)],
                            vector::zero
                        )
                    );
                    dts.append(labelledPoint(bf.size(), vector::zero));

                    DynList<labelledPoint> lpf;
                    forAll(bf, pI)
                    {
                        const labelledPoint lp
                        (
                            globalLabel[bp[bf[bpI]]],
                            points[bf[pI]]
                        );

                        dts.append(lp);
                        lpf.append(lp);
                    }

                    faceMap_[bpI].append(lpf);
                    faceAtProc_[bpI].append(Pstream::myProcNo());
                    facePatches_[bpI].append(facePatch[pointFaces(bpI, i)]);
                }
            }
        }
    }

    std::map<label, List<labelledPoint> > receivedDataMap;
    help::exchangeMap(exchangeData, receivedDataMap);

    for
    (
        std::map<label, List<labelledPoint> >::const_iterator it=receivedDataMap.begin();
        it!=receivedDataMap.end();
        ++it
    )
    {
        const label procI = it->first;
        const List<labelledPoint>& receivedData = it->second;

        for(label i=0;i<receivedData.size();)
        {
            const label bpI = globalToLocal[receivedData[i++].pointLabel()];

            const label nFaces = receivedData[i++].pointLabel();

            for(label fI=0;fI<nFaces;++fI)
            {
                const label patchI = receivedData[i++].pointLabel();
                const label size = receivedData[i++].pointLabel();

                DynList<labelledPoint> lpf(size);
                forAll(lpf, pI)
                    lpf[pI] = receivedData[i++];

                faceMap_[bpI].append(lpf);
                faceAtProc_[bpI].append(procI);
                facePatches_[bpI].append(patchI);
            }
        }
    }

    //- srot the faces in the counter-clockwise order
    for
    (
        std::map<label, DynList<label> >::iterator it=faceAtProc_.begin();
        it!=faceAtProc_.end();
        ++it
    )
    {
        const label bpI = it->first;

        DynList<label>& faceProc = it->second;
        DynList<label>& fPatches = facePatches_[bpI];
        DynList<DynList<labelledPoint, 6> >& pFaces = faceMap_[bpI];

        forAll(pFaces, i)
        {
            const DynList<labelledPoint, 6>& lpf = pFaces[i];

            label pos(-1);
            forAll(lpf, pI)
                if( lpf[i].pointLabel() == globalLabel[bpI] )
                {
                    pos = pI;
                    break;
                }

            for(label j=i+2;j<pFaces.size();++j)
            {
                const DynList<labelledPoint, 6>& olpf = pFaces[j];

                if( olpf.contains(lpf[lpf.rcIndex(pos)]) )
                {
                    label helper = faceProc[j];
                    faceProc[j] = faceProc[i+1];
                    faceProc[i+1] = helper;

                    helper = fPatches[j];
                    fPatches[j] = fPatches[i+1];
                    fPatches[i+1] = helper;

                    DynList<labelledPoint, 6> lpfHelper;
                    lpfHelper = pFaces[j];
                    pFaces[j] = pFaces[i+1];
                    pFaces[i+1] = lpfHelper;
                }
            }
        }
    }
}

void edgeExtractor::cornerEvaluator::sortedFacesAtPoint
(
    const label bpI,
    DynList<label>& pFaces
) const
{
    const meshSurfaceEngine& mse = extractor_.surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& pointFaces = mse.pointFaces();
    const VRWGraph& pointInFace = mse.pointInFaces();

    pFaces = pointFaces[bpI];

    forAll(pFaces, i)
    {
        const face& bf = bFaces[pFaces[i]];
        const label pos = pointFaces.containsAtPosition(bpI, pFaces[i]);

        const edge e = bf.faceEdge(bf.rcIndex(pointInFace(bpI, pos)));

        for(label j=i+2;j<pFaces.size();++j)
        {
            if( bFaces[pFaces[j]].which(e.start()) >= 0 )
            {
                label bfJ = pFaces[i+1];
                pFaces[i+1] = pFaces[j];
                pFaces[j] = bfJ;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

edgeExtractor::cornerEvaluator::cornerEvaluator
(
    const edgeExtractor& ee,
    const meshSurfacePartitioner& mPart
)
:
    extractor_(ee),
    partitioner_(mPart),
    faceMap_(),
    facePatches_(),
    faceAtProc_()
{
    if( Pstream::parRun() )
        createParallelAddressing();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

edgeExtractor::cornerEvaluator::~cornerEvaluator()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

void edgeExtractor::cornerEvaluator::improveCorners
(
    labelList& newBoundaryPatches
)
{
    const meshOctree& mo = extractor_.meshOctree_;

    const labelHashSet& corners = partitioner_.corners();

    const meshSurfaceEngine& mse = extractor_.surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const pointFieldPMG& points = mse.points();
    const faceList::subList& bFaces = mse.boundaryFaces();

    forAllConstIter(labelHashSet, corners, it)
    {
        const point& p = points[bPoints[it.key()]];

        if( faceMap_.find(it.key()) == faceMap_.end() )
        {
            const labelList& facePatch = extractor_.facePatch_;
            Info << "Checking corner " << it.key() << endl;

            //- get faces attached to this point
            //- sorted in the counter-clocwise order
            DynList<label> pFaces;
            sortedFacesAtPoint(it.key(), pFaces);

            Info << "Faces at corner " << pFaces << endl;

            //- find feature edges attached to the point
            DynList<labelPair> edgePatches;
            DynList<label> edgeIndex;
            forAll(pFaces, i)
            {
                const label patch0 = facePatch[pFaces[i]];
                const label patch1 = facePatch[pFaces[pFaces.fcIndex(i)]];

                if( patch0 != patch1 )
                {
                    edgeIndex.append(i);
                    edgePatches.append(labelPair(patch0, patch1));
                }
            }

            Info << "Edge patches " << edgePatches << endl;
            Info << "Edge indices " << edgeIndex << endl;

            //- find the best aligned feature edge to the feature edge
            DynList<label> bestEdgeIndices(edgePatches.size());
            forAll(edgePatches, i)
            {
                const label patch0 = edgePatches[i].first();
                const label patch1 = edgePatches[i].second();
                DynList<label> patches(2);
                patches[0] = patch0;
                patches[1] = patch1;

                //- find the proection of the first point onto the feature edge
                point mps;
                scalar dSqS;
                mo.findNearestPointToPatches(mps, dSqS, p, patches);

                scalar bestAlignment(0.0);
                label bestEdgeIndex(-1);

                forAll(pFaces, pfI)
                {
                    const face& bf = bFaces[pFaces[pfI]];

                    const label pos = bf.which(bPoints[it.key()]);

                    const point& pe = points[bf[bf.rcIndex(pos)]];

                    //- calculate the vector of current edge
                    vector ev = pe - p;
                    ev /= (mag(ev) + VSMALL);

                    //- project the endpoint onto the feature edge
                    point mpe;
                    scalar dSqE;
                    mo.findNearestPointToPatches(mpe, dSqE, pe, patches);

                    vector fv = mpe - mps;
                    fv /= (mag(fv) + VSMALL);

                    //- calculate alignment metrix between the current edge
                    //- and the feature edge
                    const scalar align = 0.5 * (1.0 + (ev & fv));

                    if( align > bestAlignment )
                    {
                        bestAlignment = align;
                        bestEdgeIndex = pfI;
                    }
                }

                bestEdgeIndices[i] = bestEdgeIndex;
            }

            Info << "\nPatches of edges at corner " << it.key()
                 << " are " << edgePatches << endl;
            Info << "Orig edge indices " << edgeIndex << endl;
            Info << "New edge indices " << bestEdgeIndices << endl;

            DynList<label> newFacePatches(pFaces.size());
            forAll(bestEdgeIndices, i)
            {
                const label patch = edgePatches[i].first();
                label index = bestEdgeIndices[i];

                while( index != bestEdgeIndices[bestEdgeIndices.fcIndex(i)] )
                {
                    newFacePatches[index] = patch;
                    index = (index + 1) % pFaces.size();
                }
            }

            bool allPatchesRemain(true);
            forAll(edgePatches, i)
                if( !newFacePatches.contains(edgePatches[i].first()) )
                {
                    allPatchesRemain = false;
                    break;
                }

            Info << "New patches at point" << newFacePatches << endl;
            if( allPatchesRemain )
            {
                forAll(pFaces, i)
                    newBoundaryPatches[pFaces[i]] = newFacePatches[i];
            }
        }
        else
        {

        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

bool edgeExtractor::findCornerCandidates()
{
    bool changed(false);

    const triSurf& surf = meshOctree_.surface();
    const pointField& sp = surf.points();

    //- create the surface partitioner and find the corners on the surface mesh
    triSurfacePartitioner surfPartitioner(surf);

    const labelList& corners = surfPartitioner.corners();

    Map<label> surfacePointToCorner;
    forAll(corners, cornerI)
        surfacePointToCorner.insert(corners[cornerI], cornerI);

    List<labelledScalar> nearestPoint
    (
        corners.size(),
        labelledScalar(0, VGREAT)
    );

    //- calculate the search range
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const pointFieldPMG& points = mesh_.points();
    const labelList& bPoints = mse.boundaryPoints();
    const labelList& bp = mse.bp();
    const faceList::subList& bFaces = mse.boundaryFaces();

    scalarList searchRange(bPoints.size(), 0.0);
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];
        const point c = bf.centre(points);

        forAll(bf, pI)
        {
            const scalar d = 2.0 * Foam::mag(c - points[bf[pI]]);
            const label bpI = bp[bf[pI]];

            searchRange[bpI] = Foam::max(d, searchRange[bpI]);
        }
    }

    if( Pstream::parRun() )
    {
        const VRWGraph& bpAtProcs = mse.beAtProcs();
        const DynList<label>& bpNeiProcs = mse.bpNeiProcs();
        const Map<label>& globalToLocal =
            mse.globalToLocalBndPointAddressing();

        std::map<label, LongList<labelledScalar> > exchangeData;
        forAll(bpNeiProcs, i)
            exchangeData.insert
            (
                std::make_pair(bpNeiProcs[i], LongList<labelledScalar>())
            );

        forAllConstIter(Map<label>, globalToLocal, iter)
        {
            const label bpI = iter();

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append
                (
                    labelledScalar(iter.key(), searchRange[bpI])
                );
            }
        }

        LongList<labelledScalar> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const label bpI = globalToLocal[receivedData[i].scalarLabel()];

            searchRange[bpI] =
                Foam::max(searchRange[bpI], receivedData[i].value());
        }
    }

    DynList<label> containedTriangles;
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40) private(containedTriangles)
    # endif
    forAll(bPoints, bpI)
    {
        const point& p = points[bPoints[bpI]];
        const scalar s = searchRange[bpI];

        boundBox bb(p-vector(s, s, s), p+vector(s, s, s));
        meshOctree_.findTrianglesInBox(bb, containedTriangles);

        //- find the nearest corner on the surface mesh
        forAll(containedTriangles, i)
        {
            const labelledTri& tri = surf[containedTriangles[i]];

            forAll(tri, pI)
            {
                const label spI = tri[pI];
                if( !surfacePointToCorner.found(spI) )
                    continue;

                const label cornerI = surfacePointToCorner[spI];

                const scalar d = Foam::magSqr(sp[spI] - points[bPoints[bpI]]);

                if( nearestPoint[cornerI].value() > d )
                {
                    //- update the nearest point to the found corner
                    # ifdef USE_OMP
                    # pragma omp critical
                    # endif
                    {
                        nearestPoint[cornerI] = labelledScalar(bpI, d);
                    }
                }
            }
        }
    }

    # ifdef DEBUGEdgeExtractor
    Info << "Nearest points to corners " << nearestPoint << endl;
    # endif

    return changed;
}

bool edgeExtractor::checkCorners()
{
    bool changed(false);

    # ifdef DEBUGEdgeExtractor
    const triSurf* surfPtr = surfaceWithPatches();
    surfPtr->writeSurface("checkCornersBefore.fms");
    deleteDemandDrivenData(surfPtr);
    # endif

    const pointFieldPMG& points = mesh_.points();
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();
    const edgeList& edges = mse.edges();
    const VRWGraph& pointFaces = mse.pointFaces();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    const triSurfacePartitioner& sPartitioner = this->partitioner();
    const List<labelHashSet>& patchPatches = sPartitioner.patchPatches();

    //- allocate a copy of boundary patches
    labelList newBoundaryPatches(facePatch_.size());

    label nCorrected;
    Map<label> otherProcNewPatch;

    label nIteration(0);

    do
    {
        nCorrected = 0;
        newBoundaryPatches = facePatch_;

        //- check whether there exist situations where a boundary face
        //- is surrounded by more faces in different patches than the
        //- faces in the current patch
        if( Pstream::parRun() )
        {
            findOtherFacePatchesParallel
            (
                otherProcNewPatch,
                &newBoundaryPatches
            );
        }

        //- update the information which edges are assigned as feature edges
        meshSurfacePartitioner mPart(mse, newBoundaryPatches);

        //- find corners in the current constelation
        const labelHashSet& corners = mPart.corners();
        const VRWGraph& pPatches = mPart.pointPatches();

        typedef std::map<label, label> mapType;
        typedef std::map<label, std::pair<point, scalar> > distMapType;

        mapType cornerIndex;
        distMapType nearestToCorner;

        # ifdef DEBUGEdgeExtractor
        Info << "Finding corners " << endl;
        # endif

        //- find nearest corners in the surface mesh
        forAllConstIter(labelHashSet, corners, it)
        {
            const label bpI = it.key();
            const DynList<label> neiPatches = pPatches[bpI];

            const point& p = points[bPoints[bpI]];
            point pMap;
            scalar dSq;
            label nsp;

            if( !meshOctree_.findNearestCorner(pMap, dSq, nsp, p, neiPatches) )
            {
                nsp = -1;
                meshOctree_.findNearestPointToPatches(pMap, dSq, p, neiPatches);
            }

            nearestToCorner[bpI] = std::make_pair(pMap, dSq);
            cornerIndex[bpI] = nsp;
        }

        std::map<label, DynList<label> > bpMappedAtSurfaceCorner;
        forAllConstIter(mapType, cornerIndex, mIt)
        {
            if( mIt->second < 0 )
            {
                //- the corner does not exist in the surface mesh
                //- this may be a situation where parts of the surface are in
                //- close proximity and are not topologically connected
                Warning << "Should not get in here. Not implemented" << endl;
            }
            else
            {
                bpMappedAtSurfaceCorner[mIt->second].append(mIt->first);
            }
        }

        # ifdef DEBUGEdgeExtractor
        for
        (
            mapType::const_iterator mIt=bpMappedAtSurfaceCorner.begin();
            mIt!=bpMappedAtSurfaceCorner.end();
            ++mIt
        )
            if( mIt->second.size() > 1 )
                Info << "Surface corner " << mIt->first
                     << " mapped mesh points " << mIt->second << endl;
        # endif

        //- for all edge nodes find their nearest counterparts in the surface
        # ifdef DEBUGEdgeExtractor
        Info << "Finding nearest edges" << endl;
        # endif

        const labelHashSet& featureEdge = mPart.featureEdges();

        distMapType nearestToEdgePoint;
        mapType edgePointIndex;
        std::map<std::pair<label, label>, labelLongList> edgesSharedByPatches;

        forAllConstIter(labelHashSet, featureEdge, it)
        {
            const label beI = it.key();

            //- get the patches bounded by this feature edge
            DynList<label> patches(2);

            if( edgeFaces.sizeOfRow(beI) == 2 )
            {
                patches[0] = newBoundaryPatches[edgeFaces(beI, 0)];
                patches[1] = newBoundaryPatches[edgeFaces(beI, 1)];
            }
            else if( edgeFaces.sizeOfRow(beI) == 1 )
            {
                patches[0] = newBoundaryPatches[edgeFaces(beI, 0)];
                patches[1] = otherProcNewPatch[beI];
            }

            //- store the edge into the right group
            std::pair<label, label> patchPair
            (
                Foam::min(patches[0], patches[1]),
                Foam::max(patches[0], patches[1])
            );
            edgesSharedByPatches[patchPair].append(beI);

            //- check if some points have already been checked
            const edge& e = edges[beI];
            const label bps = bp[e.start()];
            const label bpe = bp[e.end()];

            DynList<label> checkPoints;
            if
            (
                corners.found(bps) &&
                (nearestToEdgePoint.find(bps) == nearestToEdgePoint.end())
            )
                checkPoints.append(bps);

            if
            (
                corners.found(bpe) &&
                (nearestToEdgePoint.find(bpe) == nearestToEdgePoint.end())
            )
                checkPoints.append(bpe);

            forAll(checkPoints, i)
            {
                const point& p = points[bPoints[checkPoints[i]]];
                point pMap;
                scalar dSq;
                label nse;

                if( !meshOctree_.findNearestEdgePoint(pMap, dSq, nse, p, patches) )
                {
                    nse = -1;
                    meshOctree_.findNearestPointToPatches(pMap, dSq, p, patches);
                }

                nearestToEdgePoint[checkPoints[i]] = std::make_pair(pMap, dSq);
                edgePointIndex[checkPoints[i]] = nse;
            }
        }

        labelHashSet invalidFeatureEdges;
        for
        (
            std::map<std::pair<label, label>, labelLongList>::const_iterator mIt=edgesSharedByPatches.begin();
            mIt!=edgesSharedByPatches.end();
            ++mIt
        )
        {
            const bool validConnection =
                patchPatches[mIt->first.first].found(mIt->first.second);

            # ifdef DEBUGEdgeExtractor
            Info << "Patches " << mIt->first.first
                 << " and " << mIt->first.second
                 << " share " << mIt->second
                 << " edges " << validConnection << endl;
            # endif

            if( !validConnection )
            {
                //- find surface facets in the vicinity of the edge
                //- and check whether there exist various disconnected surface
                //- parts in the vicinity of this group of edge

                const DynList<label>& invalidEdges = mIt->second;

                forAll(invalidEdges, i)
                    invalidFeatureEdges.insert(invalidEdges[i]);
            }
        }

        # ifdef DEBUGEdgeExtractor
        Info << "Invalid feature edges " << invalidFeatureEdges << endl;
        ::exit(0);
        # endif

        //- check the vicinity of the corner point and check whether
        //- it shall be replaced by some other point in the vicinity
        std::map<label, DynList<labelPair> > facesAtCornerNewPatches;
        forAllConstIter(distMapType, nearestToCorner, it)
        {
            const label bpI = it->first;

            //- find all points connected to the corner point via a face
            labelHashSet nearPoints;
            std::map<label, DynList<label> > facesContainingPoint;
            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);
                const face& bf = bFaces[bfI];

                forAll(bf, pI)
                {
                    nearPoints.insert(bf[pI]);
                    facesContainingPoint[bf[pI]].append(bfI);
                }
            }

            # ifdef DEBUGEdgeExtractor
            forAllConstIter(labelHashSet, nearPoints, iterN)
            {
                Info << "Near point " << iterN.key() << endl;
                Info << "Faces containing near point are "
                     << facesContainingPoint[iterN.key()] << endl;
            }
            # endif

            //- find the nearest point to the location where the corner
            //- shall be mapped onto the surface mesh
            label bestPoint(-1);
            scalar minDistSq(VGREAT);
            forAllConstIter(labelHashSet, nearPoints, iter)
            {
                const point& p = points[iter.key()];

                const scalar dSq = magSqr(p - nearestToCorner[bpI].first);

                if( dSq < minDistSq )
                {
                    minDistSq = dSq;
                    bestPoint = iter.key();
                }
            }

            if( bestPoint == bPoints[bpI] )
            {
                # ifdef DEBUGEdgeExtractor
                Info << "Current corner is nearest" << endl;
                # endif

                continue;
            }

            # ifdef DEBUGEdgeExtractor
            surfPtr = surfaceWithPatches(bpI);
            surfPtr->writeSurface
            (
                "corner_"+help::scalarToText(bpI)+
                "_iter_"+help::scalarToText(nIteration)+".fms"
            );
            deleteDemandDrivenData(surfPtr);
            Info << "Best candidate " << bestPoint << endl;
            # endif

            //- sort faces and edges at the corner in counter-clockwise order
            DynList<label> pFaces, pEdges;

            DynList<label> front;
            front.append(0);

            while( front.size() )
            {
                const label fI = front.removeLastElement();
                const label bfI = pointFaces(bpI, fI);

                pFaces.append(bfI);

                const face& bf = bFaces[bfI];
                const label pos = bf.which(bPoints[bpI]);

                const label beI = faceEdges(bfI, bf.rcIndex(pos));

                pEdges.append(beI);

                label nei = edgeFaces(beI, 0);
                if( nei == bfI )
                    nei = edgeFaces(beI, 1);

                if( pFaces.contains(nei) || (nei < 0) )
                    continue;

                front.append(pointFaces.containsAtPosition(bpI, nei));
            }

            # ifdef DEBUGEdgeExtractor
            Info << "Boundary point " << bpI
                 << " pFaces " << pFaces
                 << " pEdges " << pEdges << endl;
            forAll(pFaces, i)
            {
                Info << "Face " << pFaces[i] << " has nodes " << bFaces[pFaces[i]] << endl;
                Info << "Face " << pFaces[i] << " is in patch " << facePatch_[pFaces[i]] << endl;
                Info << "Edge " << pEdges[i] << " has nodes " << edges[pEdges[i]] << endl;
            }
            # endif

            //- find the best fitting edge to move towards the corner
            label bestEdge(-1);
            scalar bestAlignment(0.0);

            vector dv = it->second.first - points[bPoints[bpI]];
            dv /= (mag(dv) + VSMALL);

            if( facesContainingPoint[bestPoint].size() == 1 )
            {
                const label bfI = facesContainingPoint[bestPoint][0];

                //- only the edges of this face can be candidates
                forAll(pEdges, i)
                {
                    const label beI = pEdges[i];
                    const edge& e = edges[pEdges[i]];

                    if( !(edgeFaces(beI, 0) == bfI || edgeFaces(beI, 1) == bfI) )
                        continue;

                    vector ev
                    (
                        points[e.otherVertex(bPoints[bpI])] - points[bPoints[bpI]]
                    );
                    ev /= (mag(ev) + VSMALL);

                    const scalar metric = 0.5 * (1.0 + (dv & ev));

                    if( metric > bestAlignment )
                    {
                        bestAlignment = metric;
                        bestEdge = pEdges[i];
                    }
                }
            }
            else
            {
                //- the edge containing the best point is the candidate
                forAll(pEdges, i)
                {
                    const edge& e = edges[pEdges[i]];

                    if( e.otherVertex(bestPoint) != -1 )
                    {
                        //- select the edge which contains best fitting point
                        bestEdge = pEdges[i];
                        break;
                    }
                }
            }

            # ifdef DEBUGEdgeExtractor
            Info << "Best edge " << bestEdge
                 << " has alignment " << bestAlignment << endl;
            # endif

            //- find groups of sorted faces at this corner which are bounded by
            //- existing feature edges and the new candidate
            DynList<label> faceInGroup(pFaces.size(), -1);
            label groupI(0);
            forAll(pFaces, i)
            {
                if( faceInGroup[i] != -1 )
                    continue;

                DynList<label> front;
                front.append(i);
                faceInGroup[i] = groupI;

                while( front.size() != 0 )
                {
                    const label currIndex = front.removeLastElement();

                    const label nbeI = pEdges[currIndex];
                    const label pbeI = pEdges[pFaces.rcIndex(currIndex)];
                    if( (nbeI != bestEdge) && !featureEdge.found(nbeI) )
                    {
                        const label nfI = pFaces.fcIndex(currIndex);

                        if( faceInGroup[nfI] == -1 )
                        {
                            faceInGroup[nfI] = groupI;
                            front.append(nfI);
                        }
                    }
                    else if( (pbeI != bestEdge) && !featureEdge.found(pbeI) )
                    {
                        const label pfI = pFaces.rcIndex(currIndex);

                        if( faceInGroup[pfI] == -1 )
                        {
                            faceInGroup[pfI] = groupI;
                            front.append(pfI);
                        }
                    }
                }

                ++groupI;
            }

            # ifdef DEBUGEdgeExtractor
            Info << "faceInGroup " << faceInGroup << endl;
            # endif

            //- check which group of faces shall change patch in order to make
            //- the best fitting edge a feature edge
            DynList<labelPair> groupPairs, groupPatches;
            labelPair groupsForChanging(-1, -1);
            forAll(pEdges, i)
            {
                const label beI = pEdges[i];

                if( beI == bestEdge )
                {
                    label pos = pFaces.containsAtPosition(edgeFaces(beI, 0));
                    groupsForChanging.first() = faceInGroup[pos];
                    pos = pFaces.containsAtPosition(edgeFaces(beI, 1));
                    groupsForChanging.second() = faceInGroup[pos];
                }
                else if( featureEdge.found(beI) )
                {
                    labelPair lp, lpp;
                    label pos = pFaces.containsAtPosition(edgeFaces(beI, 0));
                    lpp.first() = facePatch_[pFaces[pos]];
                    lp.first() = faceInGroup[pos];
                    pos = pFaces.containsAtPosition(edgeFaces(beI, 1));
                    lpp.second() = facePatch_[pFaces[pos]];
                    lp.second() = faceInGroup[pos];

                    groupPairs.append(lp);
                    groupPatches.append(lpp);
                }
            }

            # ifdef DEBUGEdgeExtractor
            Info << "group pairs " << groupPairs << endl;
            Info << "Group patches " << groupPatches << endl;
            Info << "Groups for changing " << groupsForChanging << endl;
            # endif

            //- check which groups shall change their patch
            //- desired patches at the best edge are determined by finding
            //- the best alignment between the bestEdge and the feature edges
            scalar maxAlignment(0.0);
            label bestGroupsOfPatches(-1);
            forAll(groupPatches, i)
            {
                const labelPair& pp = groupPatches[i];

                const point& bes = points[edges[bestEdge].start()];
                const point& bee = points[edges[bestEdge].end()];

                DynList<label> patches(2);
                patches[0] = pp.first();
                patches[1] = pp.second();

                point mps, mpe;
                scalar dSqS, dSqE;
                label nse;

                meshOctree_.findNearestEdgePoint(mps, dSqS, nse, bes, patches);
                meshOctree_.findNearestEdgePoint(mpe, dSqE, nse, bee, patches);

                const scalar align = 0.5 * (1.0 + ((bee - bes) & (mpe - mps)));

                if( align > maxAlignment )
                {
                    maxAlignment = align;
                    bestGroupsOfPatches = i;
                }
            }

            # ifdef DEBUGEdgeExtractor
            Info << "Best groups of patches " << bestGroupsOfPatches << endl;
            # endif

            //- assign new patches
            FixedList<label, 2> newPatchForGroup(-1);
            DynList<labelPair>& newPatchForFace = facesAtCornerNewPatches[bpI];
            forAll(groupPatches, i)
            {
                if( i == bestGroupsOfPatches )
                    continue;

                const labelPair& gp = groupPairs[i];
                const labelPair& gpp = groupPatches[i];

                if
                (
                    gp.first() == groupsForChanging.first() ||
                    gp.first() == groupsForChanging.second()
                )
                {
                    const label otherPatch = gpp.second();

                    scalar Eold(0.0), Enew(0.0);
                    forAll(faceInGroup, j)
                    {
                        if( faceInGroup[j] != gp.first() )
                            continue;

                        const label bfI = pFaces[j];

                        forAllRow(faceEdges, bfI, feI)
                        {
                            const label beI = faceEdges(bfI, feI);

                            # ifdef DEBUGEdgeExtractor
                            Info << "2. beI " << beI << endl;
                            # endif

                            const point& ps = points[edges[beI].start()];
                            const point& pe = points[edges[beI].end()];
                            const scalar magE = edges[beI].mag(points) + VSMALL;

                            vector ev = pe - ps;
                            ev /= magE;

                            label nei = edgeFaces(beI, 0);
                            if( nei == bfI )
                                nei = edgeFaces(beI, 1);

                            const label posNei = pFaces.containsAtPosition(nei);
                            if
                            (
                                (posNei < 0) ||
                                (faceInGroup[posNei] != faceInGroup[j])
                            )
                            {
                                if( facePatch_[bfI] != facePatch_[nei] )
                                {
                                    DynList<label> patches(2);
                                    patches[0] = facePatch_[bfI];
                                    patches[1] = facePatch_[nei];
                                    //- calculate deformation energy
                                    //- of the old state
                                    point mps, mpe;
                                    scalar dSqS, dSqE;

                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mps,
                                        dSqS,
                                        ps,
                                        patches
                                    );
                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mpe,
                                        dSqE,
                                        pe,
                                        patches
                                    );

                                    vector fv = mpe - mps;
                                    fv /= (mag(fv) + VSMALL);

                                    scalar c = min(fv & ev, 1.0);
                                    c = max(-1.0, c);
                                    const scalar angle = acos(c);

                                    Eold +=
                                        1.0/magE * (dSqS + dSqE) + magE * angle;
                                }
                                if( otherPatch != facePatch_[nei] )
                                {
                                    //- calculate deformation energy
                                    //- of the new state
                                    DynList<label> patches(2);
                                    patches[0] = otherPatch;
                                    patches[1] = facePatch_[nei];

                                    point mps, mpe;
                                    scalar dSqS, dSqE;

                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mps,
                                        dSqS,
                                        ps,
                                        patches
                                    );
                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mpe,
                                        dSqE,
                                        pe,
                                        patches
                                    );

                                    vector fv = mpe - mps;
                                    fv /= (mag(fv) + VSMALL);

                                    scalar c = min(fv & ev, 1.0);
                                    c = max(-1.0, c);
                                    const scalar angle = acos(c);

                                    Enew +=
                                        1.0/magE * (dSqS + dSqE) + magE * angle;
                                }
                            }
                        }
                    }

                    # ifdef DEBUGEdgeExtractor
                    Info << "1. Eold " << Eold << " Enew " << Enew << endl;
                    # endif

                    if( Enew <= Eold )
                    {
                        newPatchForGroup[0] = otherPatch;

                        forAll(faceInGroup, j)
                        {
                            if( faceInGroup[j] != gp.first() )
                                continue;

                            newPatchForFace.append
                            (
                                labelPair(pFaces[j], otherPatch)
                            );
                        }
                    }
                }
                else if
                (
                    gp.second() == groupsForChanging.first() ||
                    gp.second() == groupsForChanging.second()
                )
                {
                    const label otherPatch = gpp.first();

                    scalar Eold(0.0), Enew(0.0);
                    forAll(faceInGroup, j)
                    {
                        if( faceInGroup[j] != gp.second() )
                            continue;

                        const label bfI = pFaces[j];

                        forAllRow(faceEdges, bfI, feI)
                        {
                            const label beI = faceEdges(bfI, feI);

                            # ifdef DEBUGEdgeExtractor
                            Info << "2. beI " << beI << endl;
                            # endif

                            const point& ps = points[edges[beI].start()];
                            const point& pe = points[edges[beI].end()];
                            const scalar magE = edges[beI].mag(points) + VSMALL;

                            vector ev = pe - ps;
                            ev /= magE;

                            label nei = edgeFaces(beI, 0);
                            if( nei == bfI )
                                nei = edgeFaces(beI, 1);

                            const label posNei = pFaces.containsAtPosition(nei);
                            if
                            (
                                (posNei < 0) ||
                                (faceInGroup[posNei] != faceInGroup[j])
                            )
                            {
                                if( facePatch_[bfI] != facePatch_[nei] )
                                {
                                    DynList<label> patches(2);
                                    patches[0] = facePatch_[bfI];
                                    patches[1] = facePatch_[nei];

                                    //- calculate deformation energy
                                    //- of the old state
                                    point mps, mpe;
                                    scalar dSqS, dSqE;

                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mps,
                                        dSqS,
                                        ps,
                                        patches
                                    );
                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mpe,
                                        dSqE,
                                        pe,
                                        patches
                                    );

                                    vector fv = mpe - mps;
                                    fv /= (mag(fv) + VSMALL);

                                    scalar c = min(fv & ev, 1.0);
                                    c = max(-1.0, c);
                                    const scalar angle = acos(c);

                                    Eold +=
                                        1.0/magE * (dSqS + dSqE) + magE * angle;
                                }
                                if( otherPatch != facePatch_[nei] )
                                {
                                    //- calculate deformation energy
                                    //- of the new state
                                    DynList<label> patches(2);
                                    patches[0] = otherPatch;
                                    patches[1] = facePatch_[nei];

                                    point mps, mpe;
                                    scalar dSqS, dSqE;

                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mps,
                                        dSqS,
                                        ps,
                                        patches
                                    );
                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mpe,
                                        dSqE,
                                        pe,
                                        patches
                                    );

                                    vector fv = mpe - mps;
                                    fv /= (mag(fv) + VSMALL);

                                    scalar c = min(fv & ev, 1.0);
                                    c = max(-1.0, c);
                                    const scalar angle = acos(c);

                                    Enew +=
                                        1.0/magE * (dSqS + dSqE) + magE * angle;
                                }
                            }
                        }
                    }

                    # ifdef DEBUGEdgeExtractor
                    Info << "2. Eold " << Eold << " Enew " << Enew << endl;
                    # endif

                    if( Enew <= Eold )
                    {
                        newPatchForGroup[1] = otherPatch;

                        forAll(faceInGroup, j)
                        {
                            if( faceInGroup[j] != gp.second() )
                                continue;

                            newPatchForFace.append
                            (
                                labelPair(pFaces[j], otherPatch)
                            );
                        }
                    }
                }
            }

            # ifdef DEBUGEdgeExtractor
            Info << "New patches for boundary faces "
                 << newPatchForFace << endl;
            # endif
        }

        labelHashSet changedPatch;
        typedef std::map<label, DynList<labelPair> > labelToPairMap;
        forAllConstIter(labelToPairMap, facesAtCornerNewPatches, it)
        {
            const DynList<labelPair>& lp = it->second;
            forAll(lp, i)
            {
                if( changedPatch.found(lp[i].first()) )
                    FatalError << "Face " << lp[i].first()
                               << " is already modified" << abort(FatalError);

                changedPatch.insert(lp[i].first());

                newBoundaryPatches[lp[i].first()] = lp[i].second();
                ++nCorrected;
            }
        }

        reduce(nCorrected, sumOp<label>());

        if( nCorrected )
        {
            changed = true;
            facePatch_ = newBoundaryPatches;
        }

        if( nIteration++ > 5 )
        {
            Info << "Too many iterations" << endl;
            break;
        }

        # ifdef DEBUGEdgeExtractor
        surfPtr = surfaceWithPatches();
        surfPtr->writeSurface
        (
            "checkCornersAfterIteration_"+help::scalarToText(nIteration)+".fms"
        );
        deleteDemandDrivenData(surfPtr);
        # endif

    } while( nCorrected != 0 );

    return changed;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
