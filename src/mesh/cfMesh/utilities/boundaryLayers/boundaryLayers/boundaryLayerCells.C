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

#include "boundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "helperFunctionsPar.H"
#include "demandDrivenData.H"
#include "VRWGraphList.H"

#include "labelledPair.H"
#include "HashSet.H"

#include <map>

//#define DEBUGLayer

# ifdef DEBUGLayer
#include "polyMeshGenAddressing.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayers::createLayerCells(const labelList& patchLabels)
{
    Info << "Starting creating layer cells" << endl;

    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const edgeList& edges = mse.edges();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const labelList& faceOwners = mse.faceOwners();
    const labelList& bp = mse.bp();
    const VRWGraph& pointFaces = mse.pointFaces();

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pointPatches = mPart.pointPatches();

    //- mark patches which will be extruded into layer cells
    boolList treatPatches(mesh_.boundaries().size(), false);
    forAll(patchLabels, patchI)
    {
        const label pLabel = patchLabels[patchI];
        forAll(treatPatchesWithPatch_[pLabel], i)
            treatPatches[treatPatchesWithPatch_[pLabel][i]] = true;
    }

    //- create new faces at parallel boundaries
    const Map<label>* otherProcPatchPtr(NULL);
    const Map<label>* otherFaceProcPtr(NULL);

    if( Pstream::parRun() )
    {
        createNewFacesParallel(treatPatches);

        otherProcPatchPtr = &mse.otherEdgeFacePatch();
        otherFaceProcPtr = &mse.otherEdgeFaceAtProc();
    }

    //- create lists for new boundary faces
    VRWGraph newBoundaryFaces;
    labelLongList newBoundaryOwners;
    labelLongList newBoundaryPatches;

    //- create storage for new cells
    VRWGraphList cellsToAdd;

    //- create layer cells and store boundary faces
    const label nOldCells = mesh_.cells().size();
    forAll(bFaces, bfI)
    {
        if( treatPatches[boundaryFacePatches[bfI]] )
        {
            const face& f = bFaces[bfI];
            const label pKey = patchKey_[boundaryFacePatches[bfI]];

            DynList<DynList<label> > cellFaces;

            DynList<label> newF;

            //- store the current boundary face
            newF.clear();
            newF.append(f[0]);
            for(label pI=f.size()-1;pI>0;--pI)
                newF.append(f[pI]);
            cellFaces.append(newF);

            //- create parallel face
            forAll(f, pI)
                newF[pI] = findNewNodeLabel(f[pI], pKey);

            cellFaces.append(newF);

            newBoundaryFaces.appendList(newF);
            newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
            newBoundaryPatches.append(boundaryFacePatches[bfI]);

            //- create quad faces
            newF.setSize(4);
            forAll(f, pI)
            {
                newF[0] = f[pI];
                newF[1] = f.nextLabel(pI);
                newF[2] = findNewNodeLabel(newF[1], pKey);
                newF[3] = findNewNodeLabel(f[pI], pKey);

                cellFaces.append(newF);

                //- check if the face is at the boundary
                //- of the treated partitions
                const label edgeI = faceEdges(bfI, pI);
                if( edgeFaces.sizeOfRow(edgeI) == 2 )
                {
                    label neiFace = edgeFaces(edgeI, 0);
                    if( neiFace == bfI )
                        neiFace = edgeFaces(edgeI, 1);

                    if( !treatPatches[boundaryFacePatches[neiFace]] )
                    {
                        newBoundaryFaces.appendList(newF);
                        newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                        newBoundaryPatches.append(boundaryFacePatches[neiFace]);
                    }
                }
                else if( edgeFaces.sizeOfRow(edgeI) == 1 )
                {
                    const Map<label>& otherProcPatch = *otherProcPatchPtr;
                    if( !treatPatches[otherProcPatch[edgeI]] )
                    {
                        //- face is a new boundary face
                        newBoundaryFaces.appendList(newF);
                        newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                        newBoundaryPatches.append(otherProcPatch[edgeI]);
                    }
                }
            }

            # ifdef DEBUGLayer
            Info << "Adding cell " << cellFaces << endl;
            # endif

            cellsToAdd.appendGraph(cellFaces);
        }
        else
        {
            # ifdef DEBUGLayer
            Info << "Storing original boundary face "
                << bfI << " into patch " << boundaryFacePatches[bfI] << endl;
            # endif

            newBoundaryFaces.appendList(bFaces[bfI]);
            newBoundaryOwners.append(faceOwners[bfI]);
            newBoundaryPatches.append(boundaryFacePatches[bfI]);
        }
    }

    //- data for parallel execution
    boolList procPoint;
    LongList<DynList<label, 4> > pointProcFaces;
    LongList<labelPair> faceAtPatches;
    if( Pstream::parRun() )
    {
        procPoint.setSize(nPoints_);
        procPoint = false;

        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();
        const labelList& bPoints = mse.boundaryPoints();

        for
        (
            Map<label>::const_iterator iter=globalToLocal.begin();
            iter!=globalToLocal.end();
            ++iter
        )
        {
            const label bpI = iter();
            procPoint[bPoints[bpI]] = true;
        }
    }

    //- create cells at edges
    forAll(edgeFaces, edgeI)
    {
        //- do not consider edges with no faces attached to it
        if( edgeFaces.sizeOfRow(edgeI) == 0 )
            continue;

        //- cells are generated at the processor with the lowest label
        if(
            (edgeFaces.sizeOfRow(edgeI) == 1) &&
            (otherFaceProcPtr->operator[](edgeI) < Pstream::myProcNo())
        )
            continue;

        //- check if the edge is a feature edge
        const label patchI = boundaryFacePatches[edgeFaces(edgeI, 0)];

        label patchJ;
        if( otherProcPatchPtr && otherProcPatchPtr->found(edgeI) )
        {
            patchJ = otherProcPatchPtr->operator[](edgeI);
        }
        else
        {
            patchJ = boundaryFacePatches[edgeFaces(edgeI, 1)];
        }

        if( patchI == patchJ )
            continue;

        //- check if the faces attached to the edge have different keys
        const label pKeyI = patchKey_[patchI];
        const label pKeyJ = patchKey_[patchJ];

        if( pKeyI < 0 || pKeyJ < 0 )
        {
            continue;
            FatalErrorIn
            (
                "void boundaryLayers::createLayerCells(const labelList&)"
            ) << "Patch key is negative at concave edge" << abort(FatalError);
        }

        if( pKeyI == pKeyJ )
            continue;

        const edge& e = edges[edgeI];
        if( otherVrts_.find(e.start()) == otherVrts_.end() )
            continue;
        if( otherVrts_.find(e.end()) == otherVrts_.end() )
            continue;

        //- generate faces of the bnd layer cell
        FixedList<FixedList<label, 4>, 6> cellFaces;
        createNewCellFromEdge(e, pKeyI, pKeyJ, cellFaces);

        //- store boundary faces
        newBoundaryFaces.appendList(cellFaces[1]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchJ);

        newBoundaryFaces.appendList(cellFaces[3]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchI);

        //- check if face 5 is a boundary face or at an inter-processor boundary
        const label bps = bp[e.start()];
        label unusedPatch(-1);
        forAllRow(pointPatches, bps, i)
        {
            const label ptchI = pointPatches(bps, i);

            if( ptchI == patchI )
                continue;
            if( ptchI == patchJ )
                continue;
            if( unusedPatch != -1 )
            {
                unusedPatch = -1;
                break;
            }

            unusedPatch = ptchI;
        }

        if( unusedPatch != -1 && treatedPatch_[unusedPatch] )
        {
            //- add a face in the empty patch in case of 2D layer generation
            newBoundaryFaces.appendList(cellFaces[5]);
            newBoundaryOwners.append(nOldCells+cellsToAdd.size());
            newBoundaryPatches.append(unusedPatch);
        }
        else if( Pstream::parRun() && procPoint[e.start()] )
        {
            //- add a face at inter-pocessor boundary
            pointProcFaces.append(cellFaces[5]);
            faceAtPatches.append(labelPair(patchI, patchJ));
        }

        //- check if face 4 is a boundary face or at an inter-processor boundary
        const label bpe = bp[e.end()];
        unusedPatch = -1;
        forAllRow(pointPatches, bpe, i)
        {
            const label ptchI = pointPatches(bpe, i);

            if( ptchI == patchI )
                continue;
            if( ptchI == patchJ )
                continue;
            if( unusedPatch != -1 )
            {
                unusedPatch = -1;
                break;
            }

            unusedPatch = ptchI;
        }

        if( unusedPatch != -1 && treatedPatch_[unusedPatch] )
        {
            //- add a face in the empty patch in case of 2D layer generation
            newBoundaryFaces.appendList(cellFaces[4]);
            newBoundaryOwners.append(nOldCells+cellsToAdd.size());
            newBoundaryPatches.append(unusedPatch);
        }
        else if( Pstream::parRun() && procPoint[e.end()] )
        {
            //- add a face at inter-pocessor boundary
            pointProcFaces.append(cellFaces[4]);
            faceAtPatches.append(labelPair(patchI, patchJ));
        }

        # ifdef DEBUGLayer
        Info << "Adding new cell at edge " << cellFaces << endl;
        # endif

        //- append cell to the queue
        cellsToAdd.appendGraph(cellFaces);
    }

    //- create cells for corner nodes
    typedef std::map<std::pair<label, label>, label> mPairToLabelType;
    typedef std::map<label, mPairToLabelType> mPointsType;
    typedef std::map<label, DynList<label, 3> > ppType;

    ppType nodePatches;
    labelHashSet parPoint;

    if( Pstream::parRun() )
    {
        const labelList& bPoints = mse.boundaryPoints();
        const VRWGraph& pProcs = mse.bpAtProcs();
        const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

        std::map<label, labelLongList> facesToSend;

        typedef std::map<label, DynList<DynList<label, 8>, 8> > ppfType;

        ppfType parPointFaces;
        ppType parPointPatches;

        forAllConstIter(mPointsType, otherVrts_, iter)
        {
            //- skip points on feature edges
            if( iter->second.size() == 2 )
                continue;

            const label bpI = bp[iter->first];

            if( pProcs.sizeOfRow(bpI) != 0 )
            {
                parPoint.insert(iter->first);

                //- point is at a parallel boundary
                label pMin = pProcs(bpI, 0);
                forAllRow(pProcs, bpI, i)
                {
                    const label prI = pProcs(bpI, i);

                    if( facesToSend.find(prI) == facesToSend.end() )
                        facesToSend.insert
                        (
                            std::make_pair(prI, labelLongList())
                        );

                    if( prI < pMin )
                        pMin = prI;
                }

                if( Pstream::myProcNo() == pMin )
                {
                    DynList<label, 3>& pPatches = parPointPatches[bpI];
                    pPatches.setSize(pointFaces.sizeOfRow(bpI));

                    DynList<DynList<label, 8>, 8>& pFaces = parPointFaces[bpI];
                    pFaces.setSize(pPatches.size());

                    forAllRow(pointFaces, bpI, pfI)
                    {
                        const label bfI = pointFaces(bpI, pfI);
                        const face& bf = bFaces[bfI];

                        pPatches[pfI] = boundaryFacePatches[bfI];

                        DynList<label, 8>& bfCopy = pFaces[pfI];
                        bfCopy.setSize(bf.size());
                        forAll(bf, pI)
                            bfCopy[pI] = globalPointLabel[bp[bf[pI]]];
                    }

                    continue;
                }

                labelLongList& stp = facesToSend[pMin];

                //- send the data to the processor with the lowest label
                //- data is flatenned as follows
                //- 1. the number of faces and global point label
                //- 2. number of points in the face
                //- 3. patch label
                //- 4. global labels of face points
                stp.append(globalPointLabel[bpI]);
                stp.append(pointFaces.sizeOfRow(bpI));
                forAllRow(pointFaces, bpI, pfI)
                {
                    const label bfI = pointFaces(bpI, pfI);
                    const face& bf = bFaces[bfI];

                    stp.append(bf.size());
                    stp.append(boundaryFacePatches[bfI]);
                    forAll(bf, pI)
                        stp.append(globalPointLabel[bp[bf[pI]]]);
                }
            }
        }

        //- exchange data with other processors
        labelLongList receivedData;
        help::exchangeMap(facesToSend, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label bpI = globalToLocal[receivedData[counter++]];
            const label nFaces = receivedData[counter++];
            for(label fI=0;fI<nFaces;++fI)
            {
                DynList<label, 8> f(receivedData[counter++]);
                parPointPatches[bpI].append(receivedData[counter++]);
                forAll(f, pI)
                    f[pI] = receivedData[counter++];
                parPointFaces[bpI].append(f);
            }
        }

        //- sort faces sharing corners at the parallel boundaries
        forAllIter(ppfType, parPointFaces, iter)
        {
            DynList<DynList<label, 8>, 8>& pFaces = iter->second;
            DynList<label, 3>& fPatches = parPointPatches[iter->first];
            const label gpI = globalPointLabel[iter->first];

            for(label i=0;i<pFaces.size();++i)
            {
                const DynList<label, 8>& bf = pFaces[i];
                const label pos = bf.containsAtPosition(gpI);
                const edge e(bf[pos], bf[bf.fcIndex(pos)]);

                for(label j=i+1;j<pFaces.size();++j)
                {
                    const DynList<label, 8>& obf = pFaces[j];
                    if( obf.contains(e.start()) && obf.contains(e.end()) )
                    {
                        DynList<label, 8> add;
                        add = pFaces[i+1];
                        pFaces[i+1] = pFaces[j];
                        pFaces[j] = add;

                        const label pAdd = fPatches[i+1];
                        fPatches[i+1] = fPatches[j];
                        fPatches[j] = pAdd;
                        break;
                    }
                }
            }

            DynList<label, 3> patchIDs;
            forAll(fPatches, fpI)
                patchIDs.appendIfNotIn(fPatches[fpI]);

            nodePatches.insert(std::make_pair(bPoints[iter->first], patchIDs));
        }
    }

    //- sort out point which are not at inter-processor boundaries
    forAllConstIter(mPointsType, otherVrts_, iter)
    {
        if( iter->second.size() == 2 )
            continue;

        if( parPoint.found(iter->first) )
            continue;

        const label bpI = bp[iter->first];

        //- ensure correct orientation
        DynList<label> pFaces(pointFaces.sizeOfRow(bpI));
        forAll(pFaces, fI)
            pFaces[fI] = pointFaces(bpI, fI);

        for(label i=0;i<pFaces.size();++i)
        {
            const face& bf = bFaces[pFaces[i]];
            const edge e = bf.faceEdge(bf.which(iter->first));

            for(label j=i+1;j<pFaces.size();++j)
            {
                const face& obf = bFaces[pFaces[j]];
                if(
                    (obf.which(e.start()) >= 0) &&
                    (obf.which(e.end()) >= 0)
                )
                {
                    const label add = pFaces[i+1];
                    pFaces[i+1] = pFaces[j];
                    pFaces[j] = add;
                    break;
                }
            }
        }

        DynList<label, 3> patchIDs;
        forAll(pFaces, patchI)
        {
            patchIDs.appendIfNotIn(boundaryFacePatches[pFaces[patchI]]);
        }

        nodePatches.insert(std::make_pair(iter->first, patchIDs));
    }

    //- create layer cells for corner nodes
    forAllIter(ppType, nodePatches, iter)
    {
        const DynList<label, 3>& patchIDs = iter->second;
        DynList<label, 3> pKeys;
        forAll(patchIDs, patchI)
        {
            const label pKey = patchKey_[patchIDs[patchI]];

            if( pKey < 0 )
                continue;

            pKeys.appendIfNotIn(pKey);
        }

        if( pKeys.size() != 3 )
            continue;

        # ifdef DEBUGLayer
        Pout << "Creating corner cell at point " << iter->first << endl;
        # endif

        FixedList<FixedList<label, 4>, 6> cellFaces;
        createNewCellFromNode(iter->first, pKeys, cellFaces);

        //- store boundary faces
        newBoundaryFaces.appendList(cellFaces[1]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchIDs[0]);

        newBoundaryFaces.appendList(cellFaces[3]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchIDs[1]);

        newBoundaryFaces.appendList(cellFaces[5]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchIDs[2]);

        if( Pstream::parRun() )
        {
            if( procPoint[iter->first] )
            {
                pointProcFaces.append(cellFaces[0]);
                faceAtPatches.append(labelPair(patchIDs[1], patchIDs[2]));

                pointProcFaces.append(cellFaces[2]);
                faceAtPatches.append(labelPair(patchIDs[0], patchIDs[2]));

                pointProcFaces.append(cellFaces[4]);
                faceAtPatches.append(labelPair(patchIDs[0], patchIDs[1]));
            }
        }

        # ifdef DEBUGLayer
        Info << "Adding corner cell " << cellFaces << endl;
        # endif

        //- append cell to the queue
        cellsToAdd.appendGraph(cellFaces);
    }

    if( Pstream::parRun() )
    {
        //- create faces at parallel boundaries created from
        //- points at parallel boundaries
        createNewFacesFromPointsParallel
        (
            pointProcFaces,
            faceAtPatches
        );
    }

    //- create mesh modifier
    polyMeshGenModifier meshModifier(mesh_);

    meshModifier.addCells(cellsToAdd);

    cellsToAdd.clear();
    meshModifier.reorderBoundaryFaces();
    meshModifier.replaceBoundary
    (
        patchNames_,
        newBoundaryFaces,
        newBoundaryOwners,
        newBoundaryPatches
    );

    PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();
    forAll(boundaries, patchI)
        boundaries[patchI].patchType() = patchTypes_[patchI];

    //- delete meshSurfaceEngine
    this->clearOut();

    Info << "Finished creating layer cells" << endl;
}

void boundaryLayers::createNewFacesFromPointsParallel
(
    const LongList<DynList<label, 4> >& faceCandidates,
    const LongList<labelPair>& candidatePatches
)
{
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const labelList& bp = mse.bp();
    const VRWGraph& bpAtProcs = mse.bpAtProcs();
    const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
    const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

    labelList otherFaceProc(faceCandidates.size(), -1);
    //- some faces may appear more than once
    //- such faces are ordinary internal faces
    VRWGraph pointFaceCandidates(nPoints_);
    forAll(faceCandidates, fI)
    {
        forAll(faceCandidates[fI], pI)
            pointFaceCandidates.append(faceCandidates[fI][pI], fI);
    }

    boolList duplicateFace(faceCandidates.size(), false);
    List<labelledPair> pointOfOrigin(faceCandidates.size());
    std::map<labelledPair, label> pointOfOriginToFaceLabel;
    forAll(faceCandidates, fI)
    {
        const DynList<label, 4>& f = faceCandidates[fI];

        const label pointI = f[0];

        const labelledPair lp
        (
            globalPointLabel[bp[pointI]],
            Pair<label>
            (
                patchKey_[candidatePatches[fI][0]],
                patchKey_[candidatePatches[fI][1]]
            )
        );

        if(
            pointOfOriginToFaceLabel.find(lp) != pointOfOriginToFaceLabel.end()
        )
        {
            duplicateFace[fI] = true;
            pointOfOrigin[fI] = lp;
            duplicateFace[pointOfOriginToFaceLabel[lp]] = true;
            continue;
        }

        pointOfOrigin[fI] = lp;

        pointOfOriginToFaceLabel.insert(std::make_pair(lp, fI));
    }

    //- find the processor patch for each processor boundary face
    //- the key of the algorithm is the point from which the face was created
    //- by sending the point label and the associated patches, it will be
    //- possible to find the other processor containing that face
    std::map<label, LongList<labelledPair> > exchangeData;
    const DynList<label>& neiProcs = mse.bpNeiProcs();
    forAll(neiProcs, procI)
    {
        const label neiProcI = neiProcs[procI];

        if( neiProcI == Pstream::myProcNo() )
            continue;

        if( exchangeData.find(neiProcI) == exchangeData.end() )
            exchangeData.insert
            (
                std::make_pair(neiProcI, LongList<labelledPair>())
            );
    }

    forAll(faceCandidates, fI)
    {
        if( duplicateFace[fI] )
            continue;

        const label bpI = bp[faceCandidates[fI][0]];

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProcNo = bpAtProcs(bpI, procI);
            if( neiProcNo == Pstream::myProcNo() )
                continue;

            LongList<labelledPair>& dataToSend = exchangeData[neiProcNo];
            dataToSend.append(pointOfOrigin[fI]);
        }
    }

    //- exchange the data with other processors
    std::map<label, List<labelledPair> > receivedMap;
    help::exchangeMap(exchangeData, receivedMap);
    exchangeData.clear();

    for
    (
        std::map<label, List<labelledPair> >::const_iterator
        iter=receivedMap.begin();
        iter!=receivedMap.end();
        ++iter
    )
    {
        const List<labelledPair>& receivedData = iter->second;

        forAll(receivedData, i)
        {
            const labelledPair& lpp = receivedData[i];
            const label gpI = lpp.pairLabel();
            const label pointI = bPoints[globalToLocal[gpI]];
            const labelPair& lp = lpp.pair();

            forAllRow(pointFaceCandidates, pointI, i)
            {
                const label fI = pointFaceCandidates(pointI, i);
                const DynList<label, 4>& f = faceCandidates[fI];

                const labelPair pk
                (
                    patchKey_[candidatePatches[fI][0]],
                    patchKey_[candidatePatches[fI][1]]
                );

                const labelPair rpk
                (
                    patchKey_[candidatePatches[fI][1]],
                    patchKey_[candidatePatches[fI][0]]
                );

                if(
                    (f[0] == pointI) && ((pk == lp) || (rpk == lp))
                )
                {
                    //- found the processor containing other face
                    otherFaceProc[pointOfOriginToFaceLabel[lpp]] = iter->first;
                }
            }
        }
    }
    receivedMap.clear();

    //- sort the points in ascending order
    //- this ensures the correct order of faces at the processor boundaries
    sort(pointOfOrigin);

    Map<label> otherProcToProcPatch;
    forAll(mesh_.procBoundaries(), patchI)
    {
        const processorBoundaryPatch& wp = mesh_.procBoundaries()[patchI];
        otherProcToProcPatch.insert(wp.neiProcNo(), patchI);
    }

    //- store processor faces
    VRWGraph newProcFaces;
    labelLongList newProc;

    forAll(pointOfOrigin, i)
    {
        const label fI = pointOfOriginToFaceLabel[pointOfOrigin[i]];

        if( duplicateFace[fI] || (otherFaceProc[fI] == -1) )
            continue;

        if( !otherProcToProcPatch.found(otherFaceProc[fI]) )
        {
            otherProcToProcPatch.insert
            (
                otherFaceProc[fI],
                polyMeshGenModifier(mesh_).addProcessorPatch
                (
                    otherFaceProc[fI]
                )
            );
        }

        newProcFaces.appendList(faceCandidates[fI]);
        newProc.append(otherProcToProcPatch[otherFaceProc[fI]]);
    }

    polyMeshGenModifier(mesh_).addProcessorFaces(newProcFaces, newProc);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
