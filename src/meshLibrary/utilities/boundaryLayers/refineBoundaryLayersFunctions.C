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

#include "refineBoundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "polyMeshGenAddressing.H"
#include "polyMeshGen2DEngine.H"
#include "VRWGraphList.H"

#include "labelledPair.H"
#include "labelledScalar.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGLayer

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace bndLayerOps
{

class meshBndLayerNeighbourOperator
{
    const meshSurfaceEngine& mse_;
    const label size_;

public:

    meshBndLayerNeighbourOperator(const meshSurfaceEngine& mse)
    :
        mse_(mse),
        size_(mse.boundaryFaces().size())
    {}

    label size() const
    {
        return size_;
    }

    void operator()(const label bfI, DynList<label>& neighbourFaces) const
    {
        neighbourFaces.clear();

        const cellListPMG& cells = mse_.cells();

        const labelList& faceOwner = mse_.faceOwners();
        const label cellI = faceOwner[bfI];
        const cell& c = cells[cellI];

        const VRWGraph& faceEdges = mse_.faceEdges();
        const VRWGraph& edgeFaces = mse_.edgeFaces();

        forAllRow(faceEdges, bfI, feI)
        {
            const label edgeI = faceEdges(bfI, feI);

            if( edgeFaces.sizeOfRow(edgeI) == 2 )
            {
                label nei = edgeFaces(edgeI, 0);

                if( nei == bfI )
                    nei = edgeFaces(edgeI, 1);

                //- faces must not be part of the same cell
                if( faceOwner[nei] == cellI )
                    continue;

                //- owner cell of the other face must
                //- have cellI as its neighbour
                const cell& neiC = cells[faceOwner[nei]];
                bool sharedFace(false);
                forAll(c, fI)
                {
                    forAll(neiC, fJ)
                    {
                        if( c[fI] == neiC[fJ] )
                        {
                            sharedFace = true;
                            break;
                        }
                    }

                    if( sharedFace )
                        break;
                }

                if( sharedFace )
                    neighbourFaces.append(nei);
            }
        }
    }

    template<class labelListType>
    void collectGroups
    (
        std::map<label, DynList<label> >& neiGroups,
        const labelListType& elementInGroup,
        const DynList<label>& localGroupLabel
    ) const
    {
        const polyMeshGen& mesh = mse_.mesh();
        const faceListPMG& faces = mesh.faces();
        const cellListPMG& cells = mesh.cells();

        const edgeList& edges = mse_.edges();
        const labelList& faceOwner = mse_.faceOwners();
        const VRWGraph& edgeFaces = mse_.edgeFaces();
        const Map<label>& otherProc = mse_.otherEdgeFaceAtProc();
        const Map<label>& globalToLocal = mse_.globalToLocalBndEdgeAddressing();

        std::map<label, LongList<labelPair> > exchangeData;
        forAll(mse_.beNeiProcs(), procI)
            exchangeData[mse_.beNeiProcs()[procI]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label beI = it();

            //- combine data if the cell attached to this face has a face
            //- attached to the inter-processor boundary
            //- this must hold for boundary layer cells
            const cell& c = cells[faceOwner[edgeFaces(beI, 0)]];

            bool validCell(false);
            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                forAll(f, eI)
                {
                    const edge fe = f.faceEdge(eI);

                    if( fe == edges[beI] && mesh.faceIsInProcPatch(c[fI]) >= 0 )
                    {
                        validCell = true;
                        break;
                    }
                }

                if( validCell )
                    break;
            }

            if( !validCell )
                continue;

            const label groupI = elementInGroup[edgeFaces(beI, 0)];

            if( groupI < 0 )
                continue;

            const label lgI = localGroupLabel[groupI];
            exchangeData[otherProc[beI]].append(labelPair(it.key(), lgI));
        }

        LongList<labelPair> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const labelPair& lp = receivedData[i];

            const label beI = globalToLocal[lp.first()];

            const cell& c = cells[faceOwner[edgeFaces(beI, 0)]];

            //- combine data if the cell attached to this face has a face
            //- attached to the inter-processor boundary
            //- this must hold for boundary layer cells
            bool validCell(false);
            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                forAll(f, eI)
                {
                    const edge fe = f.faceEdge(eI);

                    if( fe == edges[beI] && mesh.faceIsInProcPatch(c[fI]) >= 0 )
                    {
                        validCell = true;
                        break;
                    }
                }

                if( validCell )
                    break;
            }

            if( !validCell )
                continue;

            const label groupI = elementInGroup[edgeFaces(beI, 0)];

            if( groupI < 0 )
                continue;

            DynList<label>& ng = neiGroups[localGroupLabel[groupI]];

            //- store the connection over the inter-processor boundary
            ng.appendIfNotIn(lp.second());
        }
    }
};

class meshBndLayerSelectorOperator
{
    const meshSurfaceEngine& mse_;

public:

    meshBndLayerSelectorOperator(const meshSurfaceEngine& mse)
    :
        mse_(mse)
    {}

    bool operator()(const label bfI) const
    {
        const labelList& faceOwner = mse_.faceOwners();
        const polyMeshGen& mesh = mse_.mesh();
        const faceListPMG& faces = mesh.faces();

        const cell& c = mesh.cells()[faceOwner[bfI]];
        const PtrList<boundaryPatch>& boundaries = mesh.boundaries();
        const label start = boundaries[0].patchStart();

        label nBndFaces(0), baseFace(-1), otherBase(-1), nQuads(0);
        forAll(c, fI)
        {
            if( faces[c[fI]].size() == 4 )
                ++nQuads;

            if( (c[fI] - start) == bfI )
            {
                baseFace = fI;
                ++nBndFaces;
            }
        }

        if( nQuads == 6 )
        {
            //- cell is a hex
            return true;
        }

        if( (nQuads + 2) != c.size() )
            return false;

        if( nBndFaces != 1 )
            return false;

        label nQuadsAttachedToBaseFace(0);
        forAll(c, fI)
        {
            if( fI == baseFace )
                continue;

            const bool sEdge =
                help::shareAnEdge(faces[c[baseFace]], faces[c[fI]]);

            if( (faces[c[fI]].size() == 4) && sEdge )
            {
                ++nQuadsAttachedToBaseFace;
            }
            else if( !sEdge )
            {
                if( otherBase != -1 )
                    return false;

                otherBase = fI;
            }
        }

        if(
            (nQuads == 6) ||
            (
                ((nQuadsAttachedToBaseFace + 2) == c.size()) &&
                otherBase != -1 &&
                !help::shareAnEdge(faces[c[baseFace]], faces[c[otherBase]])
            )
        )
            return true;

        return false;
    }
};

} // End namespace bndLayerOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void refineBoundaryLayers::analyseLayers()
{
    Info << "Analysing mesh for bnd layer existence" << endl;

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& facePatch = mse.boundaryFacePatches();
    mse.faceOwners();
    mse.faceEdges();
    mse.edgeFaces();
    mse.edges();
    mse.boundaryPointEdges();

    //- find layers in patch
    labelLongList bndFaceInLayer;
    const label nGroups =
        help::groupMarking
        (
            bndFaceInLayer,
            bndLayerOps::meshBndLayerNeighbourOperator(mse),
            bndLayerOps::meshBndLayerSelectorOperator(mse)
        );

    # ifdef DEBUGLayer
    Info << "Number of independent layers in the mesh is " << nGroups << endl;
    # endif

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    //- create patch name to index addressing
    std::map<word, label> patchNameToIndex;
    forAll(boundaries, patchI)
        patchNameToIndex[boundaries[patchI].patchName()] = patchI;

    //- check layer labels over a patch
    List<DynList<label> > groupsAtPatch(boundaries.size());
    forAll(facePatch, bfI)
    {
        if( bndFaceInLayer[bfI] < 0 )
            continue;

        groupsAtPatch[facePatch[bfI]].appendIfNotIn(bndFaceInLayer[bfI]);
    }

    //- set the information which patches have an extruded layer
    labelList groupIDs(nGroups, -1);

    layerAtPatch_.setSize(boundaries.size());
    layerAtPatch_ = -1;

    label nValidLayers(0);
    forAll(groupsAtPatch, patchI)
    {
        if( groupsAtPatch[patchI].size() == 1 )
        {
            const label groupI = groupsAtPatch[patchI][0];

            if( groupIDs[groupI] == -1 )
                groupIDs[groupI] = nValidLayers++;

            layerAtPatch_[patchI] = groupIDs[groupI];
        }
    }

    # ifdef DEBUGLayer
    Info << "Layer at patch " << layerAtPatch_ << endl;
    # endif

    //- set the information which patches are a single boundary layer face
    patchesInLayer_.setSize(nValidLayers);
    forAll(layerAtPatch_, patchI)
    {
        if( layerAtPatch_[patchI] < 0 )
            continue;

        patchesInLayer_[layerAtPatch_[patchI]].append
        (
            boundaries[patchI].patchName()
        );
    }

    # ifdef DEBUGLayer
    Info << "Patches in layer " << patchesInLayer_ << endl;

    //- write layers to a subset
    std::map<label, label> layerId;
    for(label i=0;i<nValidLayers;++i)
        layerId[i] = mesh_.addFaceSubset("layerFaces_"+help::scalarToText(i));

    forAll(layerAtPatch_, i)
    {
        if( layerAtPatch_[i] < 0 )
            continue;

        const label start = boundaries[i].patchStart();
        const label end = start + boundaries[i].patchSize();

        for(label faceI=start;faceI<end;++faceI)
            mesh_.addFaceToSubset(layerId[layerAtPatch_[i]], faceI);
    }
    mesh_.write();
    # endif

    //- set the number of boundary layers for each patch
    labelList nLayersAtPatch(layerAtPatch_.size(), -1);
    boolList protectedValue(layerAtPatch_.size(), false);

    forAll(patchesInLayer_, layerI)
    {
        const DynList<word>& layerPatches = patchesInLayer_[layerI];

        label maxNumLayers(1);
        bool hasLocalValue(false);

        //- find the maximum requested number of layers over the layer
        forAll(layerPatches, lpI)
        {
            const word pName = layerPatches[lpI];

            std::map<word, label>::const_iterator it =
                numLayersForPatch_.find(pName);

            if( it != numLayersForPatch_.end() )
            {
                //- check if the layer is interrupted at this patch
                if(
                    discontinuousLayersForPatch_.find(pName) !=
                    discontinuousLayersForPatch_.end()
                )
                {
                    //- set the number of layers and lock this location
                    nLayersAtPatch[patchNameToIndex[pName]] = it->second;
                    protectedValue[patchNameToIndex[pName]] = true;
                    hasLocalValue = true;
                }
                else
                {
                    //- take the maximum number of layers
                    maxNumLayers = Foam::max(maxNumLayers, it->second);
                    hasLocalValue = true;
                }
            }
        }

        //- apply the global value if no local values exist
        if( !hasLocalValue )
            maxNumLayers = globalNumLayers_;

        //- apply the maximum number of ayer of all unprotected patches
        forAll(layerPatches, lpI)
        {
            const label ptchI = patchNameToIndex[layerPatches[lpI]];

            if( !protectedValue[ptchI] )
                nLayersAtPatch[ptchI] = maxNumLayers;
        }
    }

    if( is2DMesh_ )
    {
        polyMeshGen2DEngine mesh2DEngine(mesh_);
        const boolList& zMinPoint = mesh2DEngine.zMinPoints();
        const boolList& zMaxPoint = mesh2DEngine.zMaxPoints();

        const faceList::subList& bFaces = mse.boundaryFaces();

        boolList allZMax(mesh_.boundaries().size(), true);
        boolList allZMin(mesh_.boundaries().size(), true);

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            forAll(bf, pI)
            {
                if( !zMinPoint[bf[pI]] )
                    allZMin[facePatch[bfI]] = false;
                if( !zMaxPoint[bf[pI]] )
                    allZMax[facePatch[bfI]] = false;
            }
        }

        //- mark empty patches as already used
        forAll(allZMin, patchI)
        {
            if( allZMin[patchI] ^ allZMax[patchI] )
            {
                nLayersAtPatch[patchI] = -1;
                layerAtPatch_[patchI] = -1;
            }
        }
    }

    # ifdef DEBUGLayer
    Pout << "nLayersAtPatch " << nLayersAtPatch << endl;
    # endif

    //- set the number of boundary layers which shall be generated above
    //- each boundary face
    nLayersAtBndFace_.setSize(facePatch.size());
    nLayersAtBndFace_ = globalNumLayers_;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(nLayersAtBndFace_, bfI)
    {
        const label patchI = facePatch[bfI];

        if( nLayersAtPatch[patchI] < 0 )
        {
            nLayersAtBndFace_[bfI] = 1;
        }
        else
        {
            nLayersAtBndFace_[bfI] = nLayersAtPatch[patchI];
        }
    }

    # ifdef DEBUGLayer
    forAll(nLayersAtBndFace_, bfI)
    Pout << "Boundary face " << bfI << " in patch "
        << facePatch[bfI] << " num layers " << nLayersAtBndFace_[bfI] << endl;
    //::exit(1);
    # endif
}

void refineBoundaryLayers::calculateAddressing
(
    const label bfI,
    label& baseFace,
    DynList<edge, 48>& edges,
    DynList<DynList<label, 2>, 48>& edgeFaces,
    DynList<DynList<label, 10>, 24>& faceEdges
) const
{
    const label nInternalFaces = mesh_.boundaries()[0].patchStart();

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& faceOwner = mse.faceOwners();

    const faceListPMG& faces = mesh_.faces();
    const cell& c = mesh_.cells()[faceOwner[bfI]];

    faceEdges.setSize(c.size());
    baseFace = -1;
    forAll(c, fI)
    {
        if( c[fI] - nInternalFaces == bfI )
        {
            baseFace = fI;
        }

        const face& f = faces[c[fI]];
        faceEdges[fI].setSize(f.size());

        forAll(f, eI)
        {
            const edge e = f.faceEdge(eI);

            label pos = edges.containsAtPosition(e);

            if( pos < 0 )
            {
                pos = edges.size();
                edges.append(e);
                edgeFaces.setSize(pos+1);
            }

            edgeFaces[pos].append(fI);
            faceEdges[fI][eI] = pos;
        }
    }
}

bool refineBoundaryLayers::findHairsForFace
(
    const label bfI,
    DynList<edge>& hairEdges
) const
{
    const label nInternalFaces = mesh_.boundaries()[0].patchStart();

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& faceOwner = mse.faceOwners();

    const faceListPMG& faces = mesh_.faces();
    const cell& c = mesh_.cells()[faceOwner[bfI]];

    //- check cell topology
    DynList<edge, 48> edges;
    DynList<DynList<label, 2>, 48> edgeFaces;
    DynList<DynList<label, 10>, 24> faceEdges;
    faceEdges.setSize(c.size());
    label baseFace(-1);
    forAll(c, fI)
    {
        if( c[fI] - nInternalFaces == bfI )
        {
            baseFace = fI;
        }

        const face& f = faces[c[fI]];
        faceEdges[fI].setSize(f.size());

        forAll(f, eI)
        {
            const edge e = f.faceEdge(eI);

            label pos = edges.containsAtPosition(e);

            if( pos < 0 )
            {
                pos = edges.size();
                edges.append(e);
                edgeFaces.setSize(pos+1);
            }

            edgeFaces[pos].append(fI);
            faceEdges[fI][eI] = pos;
        }
    }

    if( (baseFace < 0) || ((c.size() - faces[c[baseFace]].size()) != 2) )
        return false;

    //- check if all faces attached to the base face are quads
    bool isPrism(true);

    const face& bf = faces[c[baseFace]];
    forAll(bf, pI)
    {
        const label nextEdge = faceEdges[baseFace][pI];
        const label prevEdge = faceEdges[baseFace][(pI+bf.size()-1)%bf.size()];

        if( edgeFaces[nextEdge].size() != 2 || edgeFaces[prevEdge].size() != 2 )
        {
            isPrism = false;
            break;
        }

        //- find the face attached to the edge after the current point
        label otherNextFace = edgeFaces[nextEdge][0];
        if( otherNextFace == baseFace )
            otherNextFace = edgeFaces[nextEdge][1];

        //- find the face attached to the edge before the current point
        label otherPrevFace = edgeFaces[prevEdge][0];
        if( otherPrevFace == baseFace )
            otherPrevFace = edgeFaces[prevEdge][1];

        label commonEdge;
        for(commonEdge=0;commonEdge<edges.size();++commonEdge)
            if(
                edgeFaces[commonEdge].contains(otherNextFace) &&
                edgeFaces[commonEdge].contains(otherPrevFace)
            )
                break;

        if( commonEdge == edges.size() )
        {
            isPrism = false;
            break;
        }

        //- there exists a common edge which shall be used as a hair
        if( edges[commonEdge].start() == bf[pI] )
        {
            hairEdges.append(edges[commonEdge]);
        }
        else
        {
            hairEdges.append(edges[commonEdge].reverseEdge());
        }
    }

    return isPrism;
}

bool refineBoundaryLayers::findSplitEdges()
{
    bool validLayer(true);

    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& facePatch = mse.boundaryFacePatches();
    mse.faceOwners();
    const VRWGraph& pFaces = mse.pointFaces();
    const labelList& bp = mse.bp();

    # ifdef USE_OMP
    # pragma omp parallel if( bFaces.size() > 1000 )
    # endif
    {
        edgeLongList localEdges;

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(bFaces, bfI)
        {
            if( layerAtPatch_[facePatch[bfI]] < 0 )
                continue;

            //- find hair edges for this face
            DynList<edge> hairEdges;
            if( !findHairsForFace(bfI, hairEdges) )
            {
                validLayer = false;
                continue;
            }

            const face& bf = bFaces[bfI];
            forAll(bf, pI)
            {
                //- check if every the hair shall be store or not
                //- only a hair edge from a face with the smallest label
                //- out of all faces at a points is stored
                const label bpI = bp[bf[pI]];

                bool store(true);
                forAllRow(pFaces, bpI, pfI)
                {
                    const face& obf = bFaces[pFaces(bpI, pfI)];

                    if(
                        (obf.which(hairEdges[pI].end()) < 0) &&
                        (pFaces(bpI, pfI) < bfI)
                    )
                    {
                        store = false;
                        break;
                    }
                }

                if( store )
                {
                    //- hair edge shall be stored
                    localEdges.append(hairEdges[pI]);
                }
            }
        }

        # ifdef USE_OMP
        //- find the starting element for this thread
        label startEl;
        # pragma omp critical
        {
            startEl = splitEdges_.size();

            splitEdges_.setSize(startEl+localEdges.size());
        }

        //- copy the local data to splitEdges_
        forAll(localEdges, i)
            splitEdges_[startEl++] = localEdges[i];
        # else
        //- just transfer the data to splitEdges_
        splitEdges_.transfer(localEdges);
        # endif
    }

    //- create point to split edges addressing
    splitEdgesAtPoint_.reverseAddressing(splitEdges_);

    reduce(validLayer, minOp<bool>());

    # ifdef DEBUGLayer
    for(label procI=0;procI<Pstream::nProcs();++procI)
    {
        if( procI == Pstream::myProcNo() )
        {
            Pout << "Generated split edges " << splitEdges_ << endl;
        }

        returnReduce(1, sumOp<label>());
    }
    # endif

    return validLayer;
}

void refineBoundaryLayers::generateNewVertices()
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    pointFieldPMG& points = mesh_.points();

    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& pointFaces = mse.pointFaces();
    const labelList& facePatch = mse.boundaryFacePatches();
    const labelList& bp = mse.bp();

    //- allocate the data from storing parameters applying to a split edge
    LongList<scalar> firstLayerThickness(splitEdges_.size());
    LongList<scalar> thicknessRatio(splitEdges_.size());
    labelLongList nNodesAtEdge(splitEdges_.size());

    //- count the number of vertices for each split edge
    # ifdef USE_OMP
    const label nThreads = 3 * omp_get_num_procs();
    # else
    const label nThreads = 1;
    # endif

    # ifdef USE_OMP
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
//        # ifdef USE_OMP
//        const label threadI = omp_get_thread_num();
//        # else
//        const label threadI(0);
//        # endif

        //- start counting vertices at each thread
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            //- get the requested number of boundary layers
            label nLayers(1);
            scalar ratio(globalThicknessRatio_);
            scalar thickness(globalMaxThicknessFirstLayer_);
            bool overridenThickness(false);

            const label bpI = bp[e.start()];

            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);
                const label pos = help::positionOfEdgeInFace(e, bFaces[bfI]);
                if( pos >= 0 )
                    continue;

                const word& patchName =
                    boundaries[facePatch[bfI]].patchName();

                //- overrride the global value with the maximum number of layers
                //- at this edge
                nLayers = Foam::max(nLayers, nLayersAtBndFace_[bfI]);

                //- override with the maximum ratio
                const std::map<word, scalar>::const_iterator rIt =
                    thicknessRatioForPatch_.find(patchName);
                if( rIt != thicknessRatioForPatch_.end() )
                {
                    ratio = rIt->second;
                }

                //- override with the minimum thickness set for this edge
                const std::map<word, scalar>::const_iterator tIt =
                    maxThicknessForPatch_.find(patchName);
                if( tIt != maxThicknessForPatch_.end() )
                {
                    if( overridenThickness )
                    {
                        thickness = Foam::min(thickness, tIt->second);
                    }
                    else
                    {
                        thickness = tIt->second;
                        overridenThickness = true;
                    }
                }
            }

            //- store the information
            firstLayerThickness[seI] = thickness;
            thicknessRatio[seI] = ratio;
            nNodesAtEdge[seI] = nLayers + 1;
        }
    }

    if( Pstream::parRun() )
    {
        //- transfer the information over all processor for edges
        //- at inter-processor boundaries
        const labelLongList& globalEdgeLabel =
            mesh_.addressingData().globalEdgeLabel();
        const VRWGraph& edgeAtProcs = mesh_.addressingData().edgeAtProcs();
        const Map<label>& globalToLocal =
            mesh_.addressingData().globalToLocalEdgeAddressing();
        const DynList<label>& neiProcs = mesh_.addressingData().edgeNeiProcs();
        const edgeList& edges = mesh_.addressingData().edges();
        const VRWGraph& pointEdges = mesh_.addressingData().pointEdges();

        //- exchange point number of layers
        std::map<label, LongList<labelPair> > exchangeNumLayers;
        std::map<label, LongList<labelledScalar> > exchangeThickness;
        std::map<label, LongList<labelledScalar> > exchangeRatio;
        forAll(neiProcs, i)
        {
            exchangeNumLayers.insert
            (
                std::make_pair(neiProcs[i], LongList<labelPair>())
            );
            exchangeThickness.insert
            (
                std::make_pair(neiProcs[i], LongList<labelledScalar>())
            );
            exchangeRatio.insert
            (
                std::make_pair(neiProcs[i], LongList<labelledScalar>())
            );
        }

        //- exchange the number of layers
        forAll(splitEdges_, seI)
        {
            const edge& se = splitEdges_[seI];

            const label s = se.start();
            label edgeI(-1);
            forAllRow(pointEdges, s, peI)
            {
                const label eI = pointEdges(s, peI);

                if( edges[eI] == se )
                {
                    edgeI = eI;
                    break;
                }
            }

            const label geI = globalEdgeLabel[edgeI];

            if( globalToLocal.found(geI) )
            {
                forAllRow(edgeAtProcs, edgeI, i)
                {
                    const label neiProc = edgeAtProcs(edgeI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeNumLayers[neiProc].append
                    (
                        labelPair(geI, nNodesAtEdge[seI])
                    );
                    exchangeThickness[neiProc].append
                    (
                        labelledScalar(geI, firstLayerThickness[seI])
                    );
                    exchangeRatio[neiProc].append
                    (
                        labelledScalar(geI, thicknessRatio[seI])
                    );
                }
            }
        }

        //- exchange number of layers
        LongList<labelPair> receivedNumLayers;
        help::exchangeMap(exchangeNumLayers, receivedNumLayers);

        forAll(receivedNumLayers, i)
        {
            const labelPair& lp = receivedNumLayers[i];
            const label eI = globalToLocal[lp.first()];
            const edge& e = edges[eI];
            label seI(-1);
            forAllRow(splitEdgesAtPoint_, e.start(), i)
            {
                const label seJ = splitEdgesAtPoint_(e.start(), i);
                if( splitEdges_[seJ] == e )
                {
                    seI = seJ;
                    break;
                }
            }
            nNodesAtEdge[seI] = std::max(nNodesAtEdge[seI], lp.second());
        }

        //- exchange thickness ratio
        LongList<labelledScalar> receivedScalar;
        help::exchangeMap(exchangeRatio, receivedScalar);

        forAll(receivedScalar, i)
        {
            const labelledScalar& ls = receivedScalar[i];
            const label eI = globalToLocal[ls.scalarLabel()];
            const edge& e = edges[eI];
            label seI(-1);
            forAllRow(splitEdgesAtPoint_, e.start(), i)
            {
                const label seJ = splitEdgesAtPoint_(e.start(), i);
                if( splitEdges_[seJ] == e )
                {
                    seI = seJ;
                    break;
                }
            }
            thicknessRatio[seI] = std::max(thicknessRatio[seI], ls.value());
        }

        //- exchange maximum thickness of the first layer
        receivedScalar.clear();
        help::exchangeMap(exchangeThickness, receivedScalar);

        forAll(receivedScalar, i)
        {
            const labelledScalar& ls = receivedScalar[i];
            const label eI = globalToLocal[ls.scalarLabel()];
            const edge& e = edges[eI];
            label seI(-1);
            forAllRow(splitEdgesAtPoint_, e.start(), i)
            {
                const label seJ = splitEdgesAtPoint_(e.start(), i);
                if( splitEdges_[seJ] == e )
                {
                    seI = seJ;
                    break;
                }
            }
            firstLayerThickness[seI] =
                std::min(firstLayerThickness[seI], ls.value());
        }
    }

    //- calculate the number of additional vertices which will be generated
    //- on edges of the mesh
    DynList<label> numPointsAtThread;
    numPointsAtThread.setSize(nThreads);
    numPointsAtThread = 0;

    # ifdef USE_OMP
    # pragma omp parallel for num_threads(nThreads) schedule(static, 1)
    # endif
    forAll(nNodesAtEdge, seI)
    {
        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI(0);
        # endif

        numPointsAtThread[threadI] += nNodesAtEdge[seI] - 2;
    }

    //- allocate the space in a graph storing ids of points on a split edge
    newVerticesForSplitEdge_.setSizeAndRowSize(nNodesAtEdge);

    //- calculate the number of points which will be generated
    //- on split edges
    label numPoints = points.size();
    forAll(numPointsAtThread, threadI)
    {
        const label nPts = numPointsAtThread[threadI];
        numPointsAtThread[threadI] = numPoints;
        numPoints += nPts;
    }

    points.setSize(numPoints);

    # ifdef DEBUGLayer
    Info << "Generating split vertices" << endl;
    # endif

    //- generate vertices on split edges
    # ifdef USE_OMP
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI(0);
        # endif

        label& nPoints = numPointsAtThread[threadI];

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            const vector v = e.vec(points);
            const scalar magv = mag(v);

            const label nLayers = newVerticesForSplitEdge_.sizeOfRow(seI) - 1;

            scalar firstThickness = magv / nLayers;
            if( thicknessRatio[seI] > (1. + SMALL) )
            {
                firstThickness =
                    magv /
                    (
                        (1 - Foam::pow(thicknessRatio[seI], nLayers)) /
                        (1.0 - thicknessRatio[seI])
                    );

                # ifdef DEBUGLayer
                Pout << "Thread " << threadI << endl;
                Pout << "Generating vertices at split edge "
                     << " start point " << points[e.start()]
                     << " end point " << points[e.end()] << endl;
                Pout << "Edge length " << magv << endl;
                Pout << "Thickness of the first layer "
                     << firstThickness << endl;
                # endif
            }

            firstThickness =
                Foam::min
                (
                    Foam::max(firstLayerThickness[seI], SMALL),
                    firstThickness
                );

            //- generate vertices for this edge
            newVerticesForSplitEdge_(seI, 0) = e.start();

            scalar param = firstThickness;
            const vector vec = v / (magv + VSMALL);

            for(label pI=1;pI<nLayers;++pI)
            {
                //- generate the new vertex
                const point newP = points[e.start()] + param * vec;

                # ifdef DEBUGLayer
                Pout << "Split edge " << seI << " edge points " << e
                    << " start point " << points[e.start()]
                    << " end point " << points[e.end()]
                    << " param " << param
                    << " new point " << nPoints
                    << " has coordinates " << newP << endl;
                # endif

                param += firstThickness * Foam::pow(thicknessRatio[seI], pI);

                newVerticesForSplitEdge_(seI, pI) = nPoints;
                points[nPoints++] = newP;
            }

            newVerticesForSplitEdge_(seI, nLayers) = e.end();
        }
    }

    # ifdef DEBUGLayer
    for(label procI=0;procI<Pstream::nProcs();++procI)
    {
        if( procI == Pstream::myProcNo() )
        {
            forAll(splitEdges_, seI)
            {
                Pout << "\nSplit edge " << seI << " nodes " << splitEdges_[seI]
                    << " coordinates " << points[splitEdges_[seI][0]]
                    << " " << points[splitEdges_[seI][1]]
                    << " has new points "
                    << newVerticesForSplitEdge_[seI] << endl;

                forAllRow(newVerticesForSplitEdge_, seI, i)
                    Pout << "Point " << i << " on edge ha coordinates "
                         << points[newVerticesForSplitEdge_(seI, i)] << endl;
            }
        }

        returnReduce(1, sumOp<label>());
    }

    Info << "Finished generating vertices at split edges" << endl;
    //::exit(1);
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
