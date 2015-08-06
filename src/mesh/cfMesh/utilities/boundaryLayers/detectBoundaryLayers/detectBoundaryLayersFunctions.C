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

#include "detectBoundaryLayers.H"
#include "meshSurfacePartitioner.H"
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
        const Map<label>& globalToLocal =
            mse_.globalToLocalBndEdgeAddressing();

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

                    if
                    (
                        (fe == edges[beI]) &&
                        (mesh.faceIsInProcPatch(c[fI]) >= 0)
                    )
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

                    if
                    (
                        (fe == edges[beI]) &&
                        (mesh.faceIsInProcPatch(c[fI]) >= 0)
                    )
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

void detectBoundaryLayers::analyseLayers()
{
    Info << "Analysing mesh for bnd layer existence" << endl;

    const meshSurfaceEngine& mse = meshSurface_.surfaceEngine();
    const polyMeshGen& mesh = mse.mesh();
    const PtrList<boundaryPatch>& boundaries = mesh.boundaries();

    //- allocate data needed in parallel loops
    mse.faceOwners();
    mse.faceEdges();
    mse.edgeFaces();
    mse.edges();
    mse.boundaryPointEdges();

    //- find layers in patch
    nFirstLayers_ =
        help::groupMarking
        (
            layerAtBndFace_,
            bndLayerOps::meshBndLayerNeighbourOperator(mse),
            bndLayerOps::meshBndLayerSelectorOperator(mse)
        );

    # ifdef DEBUGLayer
    labelList layerSubsetId(nFirstLayers_);
    polyMeshGen& pmg = const_cast<polyMeshGen&>(mesh);
    forAll(layerSubsetId, i)
        layerSubsetId[i] = pmg.addCellSubset("bndLayer"+help::scalarToText(i));


    forAll(layerAtBndFace_, bfI)
    {
        if( layerAtBndFace_[bfI] < 0 )
            continue;

        pmg.addCellToSubset
        (
            layerSubsetId[layerAtBndFace_[bfI]],
            mse.faceOwners()[bfI]
        );
    }
    # endif

    if( is2DMesh_ )
    {
        polyMeshGen2DEngine mesh2DEngine(mse.mesh());
        const boolList& zMinPoint = mesh2DEngine.zMinPoints();
        const boolList& zMaxPoint = mesh2DEngine.zMaxPoints();

        const faceList::subList& bFaces = mse.boundaryFaces();

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            bool allZMin(true), allZMax(true);
            forAll(bf, pI)
            {
                if( !zMinPoint[bf[pI]] )
                    allZMin = false;
                if( !zMaxPoint[bf[pI]] )
                    allZMax = false;
            }

            if( allZMax ^ allZMin )
                layerAtBndFace_[bfI] = -1;
        }
    }

    //- check the number of layer found at a patch
    typedef std::map<label, DynList<label> > patchToLayerType;
    patchToLayerType patchToLayer;

    const labelList& facePatch = meshSurface_.boundaryFacePatches();

    forAll(facePatch, bfI)
    {
        patchToLayer[facePatch[bfI]].appendIfNotIn(layerAtBndFace_[bfI]);
    }

    //- all faces of a patch must be in the same layer
    layerAtPatch_.setSize(boundaries.size());
    forAll(layerAtPatch_, i)
        layerAtPatch_[i].clear();

    for
    (
        patchToLayerType::const_iterator it=patchToLayer.begin();
        it!=patchToLayer.end();
        ++it
    )
    {
        const DynList<label>& layersAtPatch = it->second;

        forAll(layersAtPatch, i)
        {
            if( layersAtPatch[i] < 0 )
            {
                layerAtPatch_[it->first].clear();
                break;
            }
            else
            {
                layerAtPatch_[it->first].append(layersAtPatch[i]);
            }
        }
    }

    //- set the layer ID to -1 for all faces where the patch is set to -1
    forAll(facePatch, bfI)
    {
        if( layerAtPatch_[facePatch[bfI]].size() == 0 )
            layerAtBndFace_[bfI] = -1;
    }

    # ifdef DEBUGLayer
    Info << "Layer at patch " << layerAtPatch_ << endl;
    forAll(layerAtBndFace_, bfI)
        if( layerAtBndFace_[bfI] < 0 )
            Info << "0.2 No layer at boundary face " << bfI << endl;
    # endif
}

bool detectBoundaryLayers::findHairsForFace
(
    const label bfI,
    DynList<edge>& hairEdges
) const
{
    const meshSurfaceEngine& mse = meshSurface_.surfaceEngine();
    const polyMeshGen& mesh = mse.mesh();

    const label nInternalFaces = mesh.boundaries()[0].patchStart();

    const labelList& faceOwner = mse.faceOwners();

    const faceListPMG& faces = mesh.faces();
    const cell& c = mesh.cells()[faceOwner[bfI]];

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

    //- check does the base face exist and is the number of faces
    //- in the cell corresponding to a prism cell
    if( (baseFace < 0) || ((c.size() - faces[c[baseFace]].size()) != 2) )
        return false;

    //- check if all faces attached to the base face are quads
    bool isPrism(true);

    const face& bf = faces[c[baseFace]];
    hairEdges.setSize(bf.size());

    forAll(bf, pI)
    {
        const label nextEdge = faceEdges[baseFace][pI];
        const label prevEdge = faceEdges[baseFace][bf.rcIndex(pI)];

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
            hairEdges[pI] = edges[commonEdge];
        }
        else
        {
            hairEdges[pI] = edges[commonEdge].reverseEdge();
        }
    }

    return isPrism;
}

void detectBoundaryLayers::generateHairEdges()
{
    hairEdges_.clear();
    hairEdgesAtBoundaryPoint_.clear();

    const meshSurfaceEngine& mse = meshSurface_.surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
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
            if( layerAtBndFace_[bfI] < 0 )
                continue;

            //- find hair edges for this face
            DynList<edge> hairEdges;
            if( !findHairsForFace(bfI, hairEdges) )
                continue;

            const face& bf = bFaces[bfI];

            forAll(bf, pI)
            {
                //- store hair edges in a list
                const edge& he = hairEdges[pI];

                if( he.start() != bf[pI] )
                    FatalErrorIn
                    (
                        "void detectBoundaryLayers::generateHairEdges()"
                    ) << "Wrong starting point" << abort(FatalError);

                localEdges.append(he);
            }
        }

        # ifdef USE_OMP
        //- find the starting element for this thread
        label startEl;
        # pragma omp critical
        {
            startEl = hairEdges_.size();

            hairEdges_.setSize(startEl+localEdges.size());
        }

        # pragma omp barrier

        //- copy the local data to splitEdges_
        forAll(localEdges, i)
            hairEdges_[startEl++] = localEdges[i];
        # else
        //- just transfer the data to splitEdges_
        hairEdges_.transfer(localEdges);
        # endif
    }

    //- filter out duplicate edges
    VRWGraph pHairEdges;
    pHairEdges.reverseAddressing(hairEdges_);

    boolList duplicateEdge(hairEdges_.size(), false);
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(pHairEdges, pointI)
    {
        forAllRow(pHairEdges, pointI, pheI)
        {
            const label heI = pHairEdges(pointI, pheI);
            const edge& he = hairEdges_[heI];

            for(label pheJ=pheI+1;pheJ<pHairEdges.sizeOfRow(pointI);++pheJ)
            {
                const label heJ = pHairEdges(pointI, pheJ);
                const edge& nhe = hairEdges_[heJ];

                if( he == nhe )
                    duplicateEdge[heJ] = true;
            }
        }
    }

    label counter(0);
    forAll(hairEdges_, heI)
    {
        if( !duplicateEdge[heI] )
        {
            if( heI > counter )
            {
                hairEdges_[counter++] = hairEdges_[heI];
            }
            else
            {
                ++counter;
            }
        }
    }

    hairEdges_.setSize(counter);

    //- create point to split edges addressing
    hairEdgesAtBoundaryPoint_.setSize(pFaces.size());

    forAll(hairEdges_, heI)
    {
        const edge& he = hairEdges_[heI];
        hairEdgesAtBoundaryPoint_.append(bp[he.start()], heI);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
