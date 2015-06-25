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

#include "meshOctreeAddressing.H"
#include "helperFunctions.H"
#include "VRWGraphSMPModifier.H"
#include "demandDrivenData.H"
#include "meshOctree.H"
#include "labelLongList.H"
#include "triSurf.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGVrt

# ifdef DEBUGVrt
#include "OFstream.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# ifdef DEBUGVrt
void writeOctreeToVTK
(
    const fileName& fName,
    const meshOctree& octree,
    const List<direction>& boxTypes,
    const direction bType
)
{
    OFstream file(fName);

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    label nBoxes(0);
    forAll(boxTypes, leafI)
        if( boxTypes[leafI] & bType )
            ++nBoxes;

    //- write points
    file << "POINTS " << 8*nBoxes << " float\n";
    forAll(boxTypes, leafI)
    {
        if( boxTypes[leafI] & bType )
        {
            FixedList<point, 8> vertices;
            octree.returnLeaf(leafI).vertices(octree.rootBox(), vertices);

            forAll(vertices, vI)
            {
                const point& p = vertices[vI];

                file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
            }
        }
    }

    //- write boxes
    file << "\nCELLS " << nBoxes
         << " " << 9 * nBoxes << nl;

    nBoxes = 0;
    forAll(boxTypes, leafI)
    {
        if( boxTypes[leafI] & bType )
        {
            const label start = 8 * nBoxes;
            file << 8 << " " << start << " " << start+1
                      << " " << start+3 << " " << start+2
                      << " " << start+4 << " " << start+5
                      << " " << start+7 << " " << start+6 << nl;

            ++nBoxes;
        }
    }

    file << nl;

    //- write cell types
    file << "CELL_TYPES " << nBoxes << nl;
    for(label i=0;i<nBoxes;++i)
        file << 12 << nl;

    file << nl;
}
# endif

void meshOctreeAddressing::createOctreePoints() const
{
    const VRWGraph& nodeLabels = this->nodeLabels();
    const boundBox& rootBox = octree_.rootBox();

    octreePointsPtr_ = new pointField(nNodes_);
    pointField& octreePoints = *octreePointsPtr_;

    const label nLeaves = nodeLabels.size();
    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 100)
    # endif
    for(label cubeI=0;cubeI<nLeaves;++cubeI)
    {
        if( nodeLabels.sizeOfRow(cubeI) == 0 )
            continue;

        FixedList<point, 8> vertices;
        const meshOctreeCubeBasic& oc = octree_.returnLeaf(cubeI);
        oc.vertices(rootBox, vertices);

        forAllRow(nodeLabels, cubeI, nI)
        {
            const label nodeI = nodeLabels(cubeI, nI);

            octreePoints[nodeI] = vertices[nI];
        }
    }
}

void meshOctreeAddressing::createNodeLabels() const
{
    const List<direction>& boxType = this->boxType();

    nodeLabelsPtr_ = new VRWGraph(octree_.numberOfLeaves());
    VRWGraph& nodeLabels = *nodeLabelsPtr_;

    //- allocate storage for node labels
    forAll(nodeLabels, leafI)
    {
        if( boxType[leafI] )
        {
            nodeLabels.setRowSize(leafI, 8);

            forAllRow(nodeLabels, leafI, i)
                nodeLabels(leafI, i) = -1;
        }
    }

    //- start creating node labels
    nNodes_ = 0;
    DynList<label> numLocalNodes;
    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        const label nThreads = omp_get_num_threads();
        const label threadI = omp_get_thread_num();
        # else
        const label nThreads = 1;
        const label threadI = 0;
        # endif

        # ifdef USE_OMP
        # pragma omp master
        # endif
        numLocalNodes.setSize(nThreads);

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        //- count the number of nodes local to each process
        label& nLocalNodes = numLocalNodes[threadI];
        nLocalNodes = 0;

        # ifdef USE_OMP
        # pragma omp for schedule(static, 100)
        # endif
        forAll(nodeLabels, leafI)
        {
            forAllRow(nodeLabels, leafI, nI)
            {
                if( nodeLabels(leafI, nI) != -1 )
                    continue;

                FixedList<label, 8> pLeaves;
                octree_.findLeavesForCubeVertex(leafI, nI, pLeaves);

                FixedList<bool, 8> validLeaf(true);
                label minLeaf(leafI);
                forAll(pLeaves, plI)
                {
                    if( pLeaves[plI] > -1 )
                    {
                        for(label i=plI+1;i<8;++i)
                            if( pLeaves[plI] == pLeaves[i] )
                            {
                                validLeaf[plI] = false;
                                validLeaf[i] = false;
                            }

                        if( !boxType[pLeaves[plI]] )
                        {
                            validLeaf[plI] = false;
                            pLeaves[plI] = -1;
                        }

                        if( validLeaf[plI] )
                            minLeaf = Foam::min(minLeaf, pLeaves[plI]);
                    }
                    else
                    {
                        validLeaf[plI] = false;
                    }
                }

                if( (minLeaf == leafI) && validLeaf[7-nI] )
                {
                    forAll(pLeaves, plI)
                        if( validLeaf[plI] )
                        {
                            //- set node labels to -2 not to repeat searches
                            nodeLabels(pLeaves[plI], (7-plI)) = -2;
                        }

                    ++nLocalNodes;
                }
            }
        }

        //- set start node for each process
        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        label startNode(0);
        for(label i=0;i<threadI;++i)
            startNode += numLocalNodes[i];

        //- start creating node labels
        # ifdef USE_OMP
        # pragma omp for schedule(static, 100)
        # endif
        forAll(nodeLabels, leafI)
        {
            forAllRow(nodeLabels, leafI, nI)
            {
                if( nodeLabels(leafI, nI) >= 0 )
                    continue;

                FixedList<label, 8> pLeaves;
                octree_.findLeavesForCubeVertex(leafI, nI, pLeaves);

                FixedList<bool, 8> validLeaf(true);
                label minLeaf(leafI);
                forAll(pLeaves, plI)
                {
                    if( pLeaves[plI] > -1 )
                    {
                        for(label i=plI+1;i<8;++i)
                            if( pLeaves[plI] == pLeaves[i] )
                            {
                                validLeaf[plI] = false;
                                validLeaf[i] = false;
                            }

                        if( !boxType[pLeaves[plI]] )
                        {
                            validLeaf[plI] = false;
                            pLeaves[plI] = -1;
                        }

                        if( validLeaf[plI] )
                            minLeaf = Foam::min(minLeaf, pLeaves[plI]);
                    }
                    else
                    {
                        validLeaf[plI] = false;
                    }
                }

                if( (minLeaf == leafI) && validLeaf[7-nI] )
                {
                    forAll(pLeaves, plI)
                        if( validLeaf[plI] )
                        {
                            //- store the vertex at the corresponding
                            //- location in the cube
                            nodeLabels(pLeaves[plI], (7-plI)) = startNode;
                        }

                    //- store vertex label
                    ++startNode;
                }
            }
        }

        //- set the number of nodes
        # ifdef USE_OMP
        # pragma omp critical
        # endif
        {
            nNodes_ = Foam::max(nNodes_, startNode);
        }
    }

    # ifdef DEBUGVrt
    List<direction> badLeaves(nodeLabels.size(), direction(0));
    forAll(nodeLabels, leafI)
        forAllRow(nodeLabels, leafI, i)
            if( nodeLabels(leafI, i) < 0 )
                badLeaves[leafI] |= 1;
    writeOctreeToVTK("badLeaves.vtk", octree_, badLeaves, 1);

    writeOctreeToVTK("meshCells.vtk", octree_, boxType, MESHCELL);
    writeOctreeToVTK("boundaryCells.vtk", octree_, boxType, BOUNDARY);

    Info << "Checking for existence of negative node labels" << endl;
    forAll(nodeLabels, leafI)
    {
        forAllRow(nodeLabels, leafI, nI)
        {
            if( nodeLabels(leafI, nI) < 0 )
            {
                FixedList<label, 8> pLeaves;
                octree_.findLeavesForCubeVertex(leafI, nI, pLeaves);

                FixedList<bool, 8> validLeaf(true);
                label minLeaf(leafI);
                forAll(pLeaves, plI)
                {
                    if( pLeaves[plI] > -1 )
                    {
                        for(label i=plI+1;i<8;++i)
                            if( pLeaves[plI] == pLeaves[i] )
                            {
                                validLeaf[plI] = false;
                                validLeaf[i] = false;
                            }

                        if( !boxType[pLeaves[plI]] )
                        {
                            validLeaf[plI] = false;
                            pLeaves[plI] = -1;
                        }

                        if( validLeaf[plI] )
                            minLeaf = Foam::min(minLeaf, pLeaves[plI]);
                    }
                    else
                    {
                        validLeaf[plI] = false;
                    }
                }

                Info << "Min leaf " << minLeaf << endl;
                Info << "Valid leaf " << validLeaf << endl;
                Info << "pLeaves " << pLeaves << endl;
                Info << "Node position " << nI << endl;

                Info << "1.Leaf " << leafI << " node labels "
                      << nodeLabels[leafI] << endl;

                forAll(validLeaf, i)
                    if( validLeaf[i] )
                        Info << "Leaf at position " << i << " has node labels "
                             << nodeLabels[pLeaves[i]]
                             << " at level "
                             << octree_.returnLeaf(pLeaves[i]).level() << endl;
            }
        }
    }
    # endif
}

void meshOctreeAddressing::createNodeLeaves() const
{
    const List<direction>& boxType = this->boxType();
    const VRWGraph& nodeLabels = this->nodeLabels();

    //- allocate nodeLeavesPtr_
    nodeLeavesPtr_ = new FRWGraph<label, 8>(nNodes_);
    FRWGraph<label, 8>& nodeLeaves = *nodeLeavesPtr_;

    boolList storedNode(nNodes_, false);
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(nodeLabels, leafI)
    {
        forAllRow(nodeLabels, leafI, nI)
        {
            const label nodeI = nodeLabels(leafI, nI);

            if( storedNode[nodeI] )
                continue;

            storedNode[nodeI] = true;

            FixedList<label, 8> pLeaves;
            octree_.findLeavesForCubeVertex(leafI, nI, pLeaves);

            forAll(pLeaves, plI)
            {
                if( pLeaves[plI] < 0 )
                    continue;

                if( !boxType[pLeaves[plI]] )
                    pLeaves[plI] = -1;
            }

            nodeLeaves.setRow(nodeI, pLeaves);
        }
    }
}

void meshOctreeAddressing::findUsedBoxes() const
{
    boxTypePtr_ = new List<direction>(octree_.numberOfLeaves(), NONE);
    List<direction>& boxType = *boxTypePtr_;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(boxType, leafI)
    {
        const meshOctreeCubeBasic& leaf = octree_.returnLeaf(leafI);

        if(
            !octree_.hasContainedTriangles(leafI) &&
            !octree_.hasContainedEdges(leafI) &&
            (leaf.cubeType() & meshOctreeCubeBasic::INSIDE)
        )
            boxType[leafI] |= MESHCELL;
    }

    if( meshDict_.found("nonManifoldMeshing") )
    {
        const bool nonManifoldMesh
        (
            readBool(meshDict_.lookup("nonManifoldMeshing"))
        );

        if( nonManifoldMesh )
        {
            # ifdef USE_OMP
            # pragma omp parallel for schedule(dynamic, 40)
            # endif
            forAll(boxType, leafI)
            {
                const meshOctreeCubeBasic& leaf = octree_.returnLeaf(leafI);
                    if( leaf.cubeType() & meshOctreeCubeBasic::UNKNOWN )
                        boxType[leafI] |= MESHCELL;
            }
        }
    }

    if( useDATABoxes_ )
    {
        Info << "Using DATA boxes" << endl;

        forAll(boxType, leafI)
        {
            if(
                octree_.hasContainedTriangles(leafI) ||
                octree_.hasContainedEdges(leafI)
            )
                boxType[leafI] |= MESHCELL;
        }

        //- do not use boxes intersecting given patches
        if( meshDict_.found("removeCellsIntersectingPatches") )
        {
            wordHashSet patchesToRemove;

            if( meshDict_.isDict("removeCellsIntersectingPatches") )
            {
                const dictionary& dict =
                    meshDict_.subDict("removeCellsIntersectingPatches");
                const wordList patchNames = dict.toc();
                forAll(patchNames, patchI)
                    patchesToRemove.insert(patchNames[patchI]);
            }
            else
            {
                wordHashSet patchesToRemoveCopy
                (
                    meshDict_.lookup("removeCellsIntersectingPatches")
                );
                patchesToRemove.transfer(patchesToRemoveCopy);
            }

            const triSurf& ts = octree_.surface();
            boolList removeFacets(ts.size(), false);

            //- remove facets in patches
            forAllConstIter(HashSet<word>, patchesToRemove, it)
            {
                const labelList matchedPatches = ts.findPatches(it.key());
                boolList activePatch(ts.patches().size(), false);
                forAll(matchedPatches, ptchI)
                    activePatch[matchedPatches[ptchI]] = true;

                forAll(ts, triI)
                {
                    if( activePatch[ts[triI].region()] )
                        removeFacets[triI] = true;
                }
            }

            //- remove facets in subsets
            forAllConstIter(HashSet<word>, patchesToRemove, it)
            {
                const label subsetID = ts.facetSubsetIndex(it.key());
                if( subsetID >= 0 )
                {
                    labelLongList facets;
                    ts.facetsInSubset(subsetID, facets);

                    forAll(facets, i)
                        removeFacets[facets[i]] = true;
                }
            }

            //- set BOUNDARY flag to boxes intersected by the given facets
            DynList<label> containedTriangles;
            forAll(boxType, leafI)
            {
                octree_.containedTriangles(leafI, containedTriangles);

                forAll(containedTriangles, i)
                {
                    if( removeFacets[containedTriangles[i]] )
                    {
                        boxType[leafI] = NONE;
                    }
                }
            }
        }
    }
    else if( meshDict_.found("keepCellsIntersectingPatches") )
    {
        wordHashSet patchesToKeep;

        if( meshDict_.isDict("keepCellsIntersectingPatches") )
        {
            const dictionary& dict =
                meshDict_.subDict("keepCellsIntersectingPatches");
            const wordList patchNames = dict.toc();

            forAll(patchNames, patchI)
                patchesToKeep.insert(patchNames[patchI]);
        }
        else
        {
            wordHashSet patchesToKeepCopy
            (
                meshDict_.lookup("keepCellsIntersectingPatches")
            );
            patchesToKeep.transfer(patchesToKeepCopy);
        }

        const triSurf& ts = octree_.surface();
        boolList keepFacets(ts.size(), false);

        //- keep facets in patches
        forAllConstIter(HashSet<word>, patchesToKeep, it)
        {
            const labelList matchedPatches = ts.findPatches(it.key());
            boolList activePatch(ts.patches().size(), false);
            forAll(matchedPatches, ptchI)
                activePatch[matchedPatches[ptchI]] = true;

            forAll(ts, triI)
            {
                if( activePatch[ts[triI].region()] )
                    keepFacets[triI] = true;
            }
        }

        //- keep facets in subsets
        forAllConstIter(wordHashSet, patchesToKeep, it)
        {
            const label subsetID = ts.facetSubsetIndex(it.key());

            if( subsetID >= 0 )
            {
                labelLongList facets;
                ts.facetsInSubset(subsetID, facets);

                forAll(facets, i)
                    keepFacets[facets[i]] = true;
            }
        }

        //- set MESHCELL flag to boxes intersected by the given facets
        DynList<label> containedTriangles;
        forAll(boxType, leafI)
        {
            octree_.containedTriangles(leafI, containedTriangles);

            forAll(containedTriangles, i)
            {
                if( keepFacets[containedTriangles[i]] )
                {
                    boxType[leafI] = MESHCELL;
                }
            }
        }
    }

    //- set BOUNDARY flag to boxes which do not have a MESHCELL flag
    DynList<label> neighs;
    # ifdef USE_OMP
    # pragma omp parallel for if( boxType.size() > 1000 ) \
    private(neighs) schedule(dynamic, 20)
    # endif
    forAll(boxType, leafI)
    {
        if( boxType[leafI] & MESHCELL )
        {
            for(label i=0;i<6;++i)
            {
                neighs.clear();
                octree_.findNeighboursInDirection(leafI, i, neighs);

                forAll(neighs, neiI)
                {
                    const label neiLabel = neighs[neiI];

                    if( neiLabel < 0 )
                        continue;

                    if( !(boxType[neiLabel] & MESHCELL) )
                        boxType[neiLabel] = BOUNDARY;
                }
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- make sure that all processors have the same information
        //- about BOUNDARY boxes
        const labelLongList& globalLeafLabel = this->globalLeafLabel();
        const VRWGraph& leafAtProcs = this->leafAtProcs();
        const Map<label>& globalLeafToLocal =
            this->globalToLocalLeafAddressing();

        std::map<label, labelLongList> exchangeData;
        forAll(octree_.neiProcs(), procI)
            exchangeData.insert
            (
                std::make_pair
                (
                    octree_.neiProcs()[procI],
                    labelLongList()
                )
            );

        forAllConstIter(Map<label>, globalLeafToLocal, iter)
        {
            const label leafI = iter();

            if( boxType[leafI] & BOUNDARY )
            {
                forAllRow(leafAtProcs, leafI, procI)
                {
                    const label neiProc = leafAtProcs(leafI, procI);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append(globalLeafLabel[leafI]);
                }
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
            boxType[globalLeafToLocal[receivedData[i]]] = BOUNDARY;
    }
}

void meshOctreeAddressing::calculateNodeType() const
{
    const FRWGraph<label, 8>& nodeLeaves = this->nodeLeaves();

    nodeTypePtr_ = new List<direction>(nNodes_, NONE);
    List<direction>& nodeType = *nodeTypePtr_;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(nodeLeaves, nodeI)
    {
        forAllRow(nodeLeaves, nodeI, nlI)
        {
            const label leafI = nodeLeaves(nodeI, nlI);

            if( leafI == -1 )
                continue;

            const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);

            if(
                (oc.cubeType() & meshOctreeCubeBasic::OUTSIDE) ||
                (oc.cubeType() & meshOctreeCubeBasic::UNKNOWN)
            )
            {
                nodeType[nodeI] |= OUTERNODE;
                break;
            }
            else if(
                !octree_.hasContainedTriangles(leafI) &&
                (oc.cubeType() &  meshOctreeCubeBasic::INSIDE)
            )
            {
                nodeType[nodeI] |= INNERNODE;
                break;
            }
        }
    }
}

void meshOctreeAddressing::createOctreeFaces() const
{
    octreeFacesPtr_ = new VRWGraph();
    octreeFacesOwnersPtr_ = new labelLongList();
    octreeFacesNeighboursPtr_ = new labelLongList();

    const VRWGraph& nodeLabels = this->nodeLabels();
    const List<direction>& boxType = this->boxType();
    this->nodeLeaves();

    label nFaces(0);
    labelList rowSizes, chunkSizes;

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        //- faces are created and stored into helper arrays, and each thread
        //- allocates its own graph for storing faces. The faces are generated
        //- by dividing the octree leaves into chunks, and distributing these
        //- chunks over the threads. There are four chunks per each thread to
        //- improve load balancing. The number of faces generated in each chunk
        //- is stored and later in used to store the faces into the octree faces
        //- graph in the correct order
        VRWGraph helperFaces;
        labelLongList helperOwner, helperNeighbour;

        # ifdef USE_OMP
        const label nThreads = omp_get_num_threads();
        const label threadI = omp_get_thread_num();
        const label nChunks = 4 * omp_get_num_threads();
        const label chunkSize = boxType.size() / nChunks + 1;
        # else
        const label nThreads(1);
        const label threadI(0);
        const label nChunks(1);
        const label chunkSize = boxType.size();
        # endif

        # ifdef USE_OMP
        # pragma omp master
        # endif
        {
            chunkSizes.setSize(nChunks);
            chunkSizes = 0;
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        for
        (
            label chunkI=threadI;
            chunkI<nChunks;
            chunkI+=nThreads
        )
        {
            const label start = chunkSize * chunkI;
            const label end = Foam::min(start+chunkSize, boxType.size());

            const label nBefore = helperFaces.size();

            for(label leafI=start;leafI<end;++leafI)
            {
                const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);

                if( boxType[leafI] & MESHCELL )
                {
                    FixedList<label, 12> edgeCentreLabel(-1);
                    for(label i=0;i<12;++i)
                        edgeCentreLabel[i] = findEdgeCentre(leafI, i);

                    for(label fI=0;fI<6;++fI)
                    {
                        DynList<label> neighbours;
                        octree_.findNeighboursInDirection
                        (
                            leafI,
                            fI,
                            neighbours
                        );

                        if( neighbours.size() != 1 )
                            continue;

                        const label nei = neighbours[0];

                        //- stop if the neighbour is on other processor
                        if( nei == meshOctreeCubeBasic::OTHERPROC )
                            continue;

                        //- create face
                        DynList<label, 8> f;
                        for(label pI=0;pI<4;++pI)
                        {
                            const label nI =
                                meshOctreeCubeCoordinates::faceNodes_[fI][pI];
                            const label feI =
                                meshOctreeCubeCoordinates::faceEdges_[fI][pI];

                            f.append(nodeLabels(leafI, nI));

                            if( edgeCentreLabel[feI] != -1 )
                                f.append(edgeCentreLabel[feI]);
                        }

                        if( nei < 0 )
                        {
                            //- face is at the boundary of the octree
                            helperFaces.appendList(f);
                            helperOwner.append(leafI);
                            helperNeighbour.append(-1);
                        }
                        else if( boxType[nei] & MESHCELL )
                        {
                            //- face is an internal face
                            if( nei > leafI )
                            {
                                helperFaces.appendList(f);
                                helperOwner.append(leafI);
                                helperNeighbour.append(nei);
                            }
                            else if
                            (
                                octree_.returnLeaf(nei).level() < oc.level()
                            )
                            {
                                //- append a reversed face
                                label i(1);
                                for(label j=f.size()-1;j>i;--j)
                                {
                                    const label add = f[j];
                                    f[j] = f[i];
                                    f[i] = add;
                                    ++i;
                                }

                                helperFaces.appendList(f);
                                helperOwner.append(nei);
                                helperNeighbour.append(leafI);
                            }
                        }
                        else if( boxType[nei] & BOUNDARY )
                        {
                            //- face is at the boundary of the mesh cells
                            helperFaces.appendList(f);
                            helperOwner.append(leafI);
                            helperNeighbour.append(nei);
                        }
                    }
                }
                else if( boxType[leafI] & BOUNDARY )
                {
                    for(label fI=0;fI<6;++fI)
                    {
                        DynList<label> neighbours;
                        octree_.findNeighboursInDirection
                        (
                            leafI,
                            fI,
                            neighbours
                        );

                        if( neighbours.size() != 1 )
                            continue;
                        const label nei = neighbours[0];
                        if( nei < 0 )
                            continue;
                        if(
                            (boxType[nei] & MESHCELL) &&
                            (octree_.returnLeaf(nei).level() < oc.level())
                        )
                        {
                            //- add a boundary face
                            const label* fNodes =
                                meshOctreeCubeCoordinates::faceNodes_[fI];
                            face cf(4);
                            for(label i=0;i<4;++i)
                            {
                                cf[i] = nodeLabels(leafI, fNodes[i]);
                            }

                            helperFaces.appendList(cf.reverseFace());
                            helperOwner.append(nei);
                            helperNeighbour.append(leafI);
                        }
                    }
                }
            }

            //- store the size of this chunk
            chunkSizes[chunkI] = helperFaces.size() - nBefore;
        }

        //- set the sizes of faces graph
        # ifdef USE_OMP
        # pragma omp critical
        # endif
        nFaces += helperFaces.size();

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp master
        # endif
        {
            rowSizes.setSize(nFaces);
            octreeFacesPtr_->setSize(nFaces);
            octreeFacesOwnersPtr_->setSize(nFaces);
            octreeFacesNeighboursPtr_->setSize(nFaces);
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        //- set the size of face graph rows and copy owners and neighbours
        for
        (
            label chunkI=threadI;
            chunkI<nChunks;
            chunkI+=nThreads
        )
        {
            label start(0), localStart(0);
            for(label i=0;i<chunkI;++i)
                start += chunkSizes[i];
            for(label i=threadI;i<chunkI;i+=nThreads)
                localStart += chunkSizes[i];

            for(label faceI=0;faceI<chunkSizes[chunkI];++faceI)
            {
                octreeFacesOwnersPtr_->operator[](start) =
                    helperOwner[localStart];
                octreeFacesNeighboursPtr_->operator[](start) =
                    helperNeighbour[localStart];
                rowSizes[start++] = helperFaces.sizeOfRow(localStart++);
            }
        }

        # ifdef USE_OMP
        # pragma omp barrier

        //- set the size of octree faces
        # pragma omp master
        # endif
        VRWGraphSMPModifier(*octreeFacesPtr_).setSizeAndRowSize(rowSizes);

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        //- copy the data into octree faces
        for
        (
            label chunkI=threadI;
            chunkI<nChunks;
            chunkI+=nThreads
        )
        {
            label start(0), localStart(0);

            for(label i=0;i<chunkI;++i)
                start += chunkSizes[i];
            for(label i=threadI;i<chunkI;i+=nThreads)
                localStart += chunkSizes[i];

            for(label faceI=0;faceI<chunkSizes[chunkI];++faceI)
            {
                for(label i=0;i<helperFaces.sizeOfRow(localStart);++i)
                    octreeFacesPtr_->operator()(start, i) =
                        helperFaces(localStart, i);

                ++start;
                ++localStart;
            }
        }
    }

    # ifdef DEBUGVrt
    List<vector> sum(octree_.numberOfLeaves(), vector::zero);
    for(label faceI=0;faceI<octreeFacesPtr_->size();++faceI)
    {
        face f(octreeFacesPtr_->sizeOfRow(faceI));
        forAll(f, pI)
            f[pI] = octreeFacesPtr_->operator()(faceI, pI);
        const vector n = f.normal(this->octreePoints());

        sum[(*octreeFacesOwnersPtr_)[faceI]] += n;
        const label nei = (*octreeFacesNeighboursPtr_)[faceI];

        if( nei < 0 )
            continue;
        sum[nei] -= n;
    }

    forAll(sum, lfI)
    {
        if
        (
            Pstream::parRun() &&
            octree_.returnLeaf(lfI).procNo() != Pstream::myProcNo()
        )
            continue;
        if( (boxType[lfI] & MESHCELL) && (mag(sum[lfI]) > SMALL) )
            Info << "Leaf " << lfI << " is not closed " << sum[lfI] << endl;
    }
    # endif
}

void meshOctreeAddressing::calculateLeafFaces() const
{
    const labelLongList& owner = octreeFaceOwner();
    const labelLongList& neighbour = octreeFaceNeighbour();

    leafFacesPtr_ = new VRWGraph(octree_.numberOfLeaves());
    VRWGraph& leafFaces = *leafFacesPtr_;

    labelList nlf(leafFaces.size(), 0);
    forAll(owner, fI)
    {
        ++nlf[owner[fI]];
        if( neighbour[fI] < 0 )
            continue;
        ++nlf[neighbour[fI]];
    }

    forAll(nlf, leafI)
        leafFaces.setRowSize(leafI, nlf[leafI]);
    nlf = 0;

    forAll(owner, fI)
    {
        leafFaces(owner[fI], nlf[owner[fI]]++) = fI;
        if( neighbour[fI] < 0 )
            continue;
        leafFaces(neighbour[fI], nlf[neighbour[fI]]++) = fI;
    }
}

void meshOctreeAddressing::calculateNodeFaces() const
{
    const VRWGraph& octreeFaces = this->octreeFaces();
    nodeFacesPtr_ = new VRWGraph(numberOfNodes());
    VRWGraph& nodeFaces = *nodeFacesPtr_;

    VRWGraphSMPModifier(nodeFaces).reverseAddressing(octreeFaces);
    nodeFaces.setSize(numberOfNodes());
}

void meshOctreeAddressing::calculateLeafLeaves() const
{
    const labelLongList& owner = octreeFaceOwner();
    const labelLongList& neighbour = octreeFaceNeighbour();

    leafLeavesPtr_ = new VRWGraph(octree_.numberOfLeaves());
    VRWGraph& leafLeaves = *leafLeavesPtr_;

    labelList nNei(leafLeaves.size(), 0);
    forAll(owner, faceI)
    {
        if( owner[faceI] < 0 )
            continue;
        if( neighbour[faceI] < 0 )
            continue;

        ++nNei[owner[faceI]];
        ++nNei[neighbour[faceI]];
    }

    forAll(nNei, leafI)
        leafLeaves.setRowSize(leafI, nNei[leafI]);

    nNei = 0;

    forAll(owner, faceI)
    {
        if( owner[faceI] < 0 )
            continue;
        if( neighbour[faceI] < 0 )
            continue;

        leafLeaves(owner[faceI], nNei[owner[faceI]]++) = neighbour[faceI];
        leafLeaves(neighbour[faceI], nNei[neighbour[faceI]]++) = owner[faceI];
    }
}

void meshOctreeAddressing::createOctreeEdges() const
{
    const VRWGraph& faces = this->octreeFaces();

    //- allocate memory for edges, face-edges addressing
    //- and node-edges addressing
    octreeEdgesPtr_ = new LongList<edge>();
    LongList<edge>& edges = *octreeEdgesPtr_;
    faceEdgesPtr_ = new VRWGraph(faces.size());
    VRWGraph& faceEdges = *faceEdgesPtr_;
    nodeEdgesPtr_ = new VRWGraph();
    VRWGraph& nodeEdges = *nodeEdgesPtr_;
    nodeEdges.setSizeAndColumnWidth(nNodes_, 6);

    forAll(faces, faceI)
    {
        faceEdges.setRowSize(faceI, faces[faceI].size());
        forAllRow(faceEdges, faceI, feI)
            faceEdges(faceI, feI) = -1;
    }

    forAll(faces, faceI)
    {
        const label nEdges = faces.sizeOfRow(faceI);

        for(label eI=0;eI<nEdges;++eI)
        {
            const edge e
            (
                faces(faceI, eI),
                faces(faceI, (eI+1)%nEdges)
            );

            label eLabel(-1);
            forAllRow(nodeEdges, e.start(), neI)
            {
                if( edges[nodeEdges(e.start(), neI)] == e )
                {
                    eLabel = nodeEdges(e.start(), neI);
                    break;
                }
            }

            if( eLabel < 0 )
            {
                //- append new edge
                faceEdges(faceI, eI) = edges.size();
                nodeEdges.append(e.start(), edges.size());
                nodeEdges.append(e.end(), edges.size());

                edges.append(e);
            }
            else
            {
                faceEdges(faceI, eI) = eLabel;
            }
        }
    }
}

void meshOctreeAddressing::calculateLeafEdges() const
{
    const VRWGraph& edgeLeaves = this->edgeLeaves();

    leafEdgesPtr_ = new VRWGraph();
    VRWGraph& leafEdges = *leafEdgesPtr_;

    VRWGraphSMPModifier(leafEdges).reverseAddressing(edgeLeaves);
    leafEdges.setSize(octree_.numberOfLeaves());
}

void meshOctreeAddressing::calculateEdgeLeaves() const
{
    const VRWGraph& edgeFaces = this->edgeFaces();
    const labelLongList& owner = this->octreeFaceOwner();
    const labelLongList& neighbour = this->octreeFaceNeighbour();

    edgeLeavesPtr_ = new VRWGraph();
    VRWGraph& edgeLeaves = *edgeLeavesPtr_;
    edgeLeaves.setSizeAndColumnWidth(edgeFaces.size(), 4);

    forAll(edgeFaces, edgeI)
    {
        forAllRow(edgeFaces, edgeI, efI)
        {
            const label fI = edgeFaces(edgeI, efI);
            const label own = owner[fI];
            const label nei = neighbour[fI];

            edgeLeaves.appendIfNotIn(edgeI, own);

            if( nei < 0 )
                continue;
            edgeLeaves.appendIfNotIn(edgeI, nei);
        }
    }
}

void meshOctreeAddressing::calculateEdgeFaces() const
{
    const VRWGraph& faceEdges = this->faceEdges();
    edgeFacesPtr_ = new VRWGraph(octreeEdges().size());
    VRWGraph& edgeFaces = *edgeFacesPtr_;

    VRWGraphSMPModifier(edgeFaces).reverseAddressing(faceEdges);
    edgeFaces.setSize(octreeEdges().size());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
