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
#include "demandDrivenData.H"
#include "helperFunctions.H"
#include "helperFunctionsPar.H"
#include "meshSurfaceCheckInvertedVertices.H"
#include "meshSurfacePartitioner.H"
#include "polyMeshGen2DEngine.H"

#include "labelledPoint.H"
#include <map>
#include <set>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGLayer

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const meshSurfaceEngine& boundaryLayers::surfaceEngine() const
{
    if( !msePtr_ )
        msePtr_ = new meshSurfaceEngine(mesh_);

    return *msePtr_;
}

const meshSurfacePartitioner& boundaryLayers::surfacePartitioner() const
{
    if( !meshPartitionerPtr_ )
        meshPartitionerPtr_ = new meshSurfacePartitioner(surfaceEngine());

    return *meshPartitionerPtr_;
}

void boundaryLayers::findPatchesToBeTreatedTogether()
{
    if( geometryAnalysed_ )
        return;

    forAll(treatPatchesWithPatch_, patchI)
        treatPatchesWithPatch_[patchI].append(patchI);

    const meshSurfaceEngine& mse = surfaceEngine();

    const pointFieldPMG& points = mesh_.points();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const edgeList& edges = mse.edges();
    const VRWGraph& eFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    //- patches must be treated together if there exist a corner where
    //- more than three patches meet
    const labelHashSet& corners = mPart.corners();
    forAllConstIter(labelHashSet, corners, it)
    {
        const label bpI = it.key();

        if( mPart.numberOfFeatureEdgesAtPoint(bpI) > 3 )
        {
            labelHashSet commonPatches;
            DynList<label> allPatches;

            forAllRow(pPatches, bpI, patchI)
            {
                const DynList<label>& tpwp =
                    treatPatchesWithPatch_[pPatches(bpI, patchI)];

                forAll(tpwp, pJ)
                {
                    if( commonPatches.found(tpwp[pJ]) )
                        continue;

                    commonPatches.insert(tpwp[pJ]);
                    allPatches.append(tpwp[pJ]);
                }
            }

            forAllRow(pPatches, bpI, patchI)
                treatPatchesWithPatch_[pPatches(bpI, patchI)] = allPatches;

            # ifdef DEBUGLayer
            Info << "Corner " << bpI << " is shared by patches "
                << pPatches[bpI] << endl;
            Info << "All patches " << allPatches << endl;
            # endif
        }
    }

    //- patches must be treated together for concave geometries
    //- edgeClassification map counts the number of convex and concave edges
    //- for a given patch. The first counts convex edges and the second counts
    //- concave ones. If the number of concave edges is of the considerable
    //- percentage, it is treated as O-topology
    meshSurfaceCheckInvertedVertices vertexCheck(mse);
    const labelHashSet& invertedVertices = vertexCheck.invertedVertices();

    std::map<std::pair<label, label>, Pair<label> > edgeClassification;
    forAll(eFaces, eI)
    {
        if( eFaces.sizeOfRow(eI) != 2 )
            continue;

        //- check if the any of the face vertices is tangled
        const edge& e = edges[eI];
        if
        (
            !is2DMesh_ &&
            (invertedVertices.found(e[0]) || invertedVertices.found(e[1]))
        )
            continue;

        const label patch0 = boundaryFacePatches[eFaces(eI, 0)];
        const label patch1 = boundaryFacePatches[eFaces(eI, 1)];
        if( patch0 != patch1 )
        {
            std::pair<label, label> pp
            (
                Foam::min(patch0, patch1),
                Foam::max(patch0, patch1)
            );
            if( edgeClassification.find(pp) == edgeClassification.end() )
                edgeClassification.insert
                (
                    std::make_pair(pp, Pair<label>(0, 0))
                );

            const face& f1 = bFaces[eFaces(eI, 0)];
            const face& f2 = bFaces[eFaces(eI, 1)];

            if
            (
                !help::isSharedEdgeConvex(points, f1, f2) ||
                (help::angleBetweenFaces(points, f1, f2) > 0.75 * M_PI)
            )
            {
                ++edgeClassification[pp].second();
            }
            else
            {
                ++edgeClassification[pp].first();
            }
        }
    }

    if( Pstream::parRun() )
    {
        const labelList& bPoints = mse.boundaryPoints();

        //- check faces over processor edges
        const labelList& globalEdgeLabel = mse.globalBoundaryEdgeLabel();
        const Map<label>& globalToLocal = mse.globalToLocalBndEdgeAddressing();

        const DynList<label>& neiProcs = mse.beNeiProcs();
        const Map<label>& otherProcPatches = mse.otherEdgeFacePatch();
        const Map<label>& otherFaceProc = mse.otherEdgeFaceAtProc();

        //- send faces sharing processor edges to other processors
        //- faces are flattened into a single contiguous array
        const labelList& bp = mse.bp();
        const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
        const Map<label>& globalPointToLocal =
            mse.globalToLocalBndPointAddressing();

        std::map<label, LongList<labelledPoint> > exchangePoints;
        forAll(neiProcs, procI)
        {
            exchangePoints.insert
            (
                std::make_pair(neiProcs[procI], LongList<labelledPoint>())
            );
        }

        //- store faces for sending
        forAllConstIter(Map<label>, otherFaceProc, it)
        {
            const label beI = it.key();

            if( eFaces.sizeOfRow(beI) == 0 )
                continue;

            const edge& e = edges[beI];

            if
            (
                !is2DMesh_ &&
                (invertedVertices.found(e[0]) || invertedVertices.found(e[1]))
            )
                continue;

            //- do not send data if the face on other processor
            //- is in the same patch
            if( otherProcPatches[beI] == boundaryFacePatches[eFaces(beI, 0)] )
                continue;

            const face& f = bFaces[eFaces(beI, 0)];

            const label neiProc = it();

            //- each face is sent as follows
            //- 1. global edge label
            //- 2. number of face nodes
            //- 3. faces nodes and vertex coordinates
            LongList<labelledPoint>& dps = exchangePoints[neiProc];
            dps.append(labelledPoint(globalEdgeLabel[beI], point()));
            dps.append(labelledPoint(f.size(), point()));
            forAll(f, pI)
            {
                dps.append
                (
                    labelledPoint
                    (
                        globalPointLabel[bp[f[pI]]],
                        points[f[pI]]
                    )
                );
            }
        }

        LongList<labelledPoint> receivedData;
        help::exchangeMap(exchangePoints, receivedData);

        //- receive faces from other processors
        Map<label> transferredPointToLocal;

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label beI =
                globalToLocal[receivedData[counter++].pointLabel()];

            DynList<label> f(receivedData[counter++].pointLabel());
            forAll(f, pI)
            {
                const labelledPoint& lp = receivedData[counter++];

                if( globalPointToLocal.found(lp.pointLabel()) )
                {
                    //- this point already exist on this processor
                    f[pI] = bPoints[globalPointToLocal[lp.pointLabel()]];
                }
                else
                {
                    //- this point does not exist on this processor
                    //- add it to the local list of points
                    //- it will be deleted when this procedure is finished
                    if( !transferredPointToLocal.found(lp.pointLabel()) )
                    {
                        //- this point has not yet been received
                        transferredPointToLocal.insert
                        (
                            lp.pointLabel(),
                            points.size()
                        );
                        mesh_.points().append(lp.coordinates());
                    }

                    f[pI] = transferredPointToLocal[lp.pointLabel()];
                }
            }

            const face& bf = bFaces[eFaces(beI, 0)];

            const label patch0 = boundaryFacePatches[eFaces(beI, 0)];
            const label patch1 = otherProcPatches[beI];

            std::pair<label, label> pp
            (
                Foam::min(patch0, patch1),
                Foam::max(patch0, patch1)
            );
            if( edgeClassification.find(pp) == edgeClassification.end() )
                edgeClassification.insert
                (
                    std::make_pair(pp, Pair<label>(0, 0))
                );

            if(
                (otherFaceProc[beI] > Pstream::myProcNo()) &&
                (
                    !help::isSharedEdgeConvex(points, bf, f) ||
                    (help::angleBetweenFaces(points, bf, f) > 0.75 * M_PI)
                )
            )
            {
                ++edgeClassification[pp].second();
            }
            else if( otherFaceProc[beI] > Pstream::myProcNo() )
            {
                ++edgeClassification[pp].first();
            }
        }

        //- set the size of points back to their original number
        mesh_.points().setSize(nPoints_);
    }

    std::map<std::pair<label, label>, Pair<label> >::const_iterator it;
    for(it=edgeClassification.begin();it!=edgeClassification.end();++it)
    {
        const std::pair<label, label>& edgePair = it->first;
        const Pair<label>& nConvexAndConcave = it->second;

        if( nConvexAndConcave.second() != 0 )
        {
            //- number of concave edges is greater than the number
            //- of the convex ones. Treat patches together.
            const label patch0 = edgePair.first;
            const label patch1 = edgePair.second;

            //- avoid adding unused patches in case of 2D meshing
            if( treatedPatch_[patch0] || treatedPatch_[patch1] )
                continue;

            treatPatchesWithPatch_[patch0].append(patch1);
            treatPatchesWithPatch_[patch1].append(patch0);
        }
    }

    if( Pstream::parRun() )
    {
        //- make sure that all processors have the same graph
        labelLongList flattenedPatches;
        forAll(treatPatchesWithPatch_, patchI)
        {
            if( treatPatchesWithPatch_[patchI].size() <= 1 )
                continue;

            flattenedPatches.append(patchI);
            flattenedPatches.append(treatPatchesWithPatch_[patchI].size());
            forAll(treatPatchesWithPatch_[patchI], itemI)
                flattenedPatches.append(treatPatchesWithPatch_[patchI][itemI]);
        }

        labelListList procPatches(Pstream::nProcs());
        procPatches[Pstream::myProcNo()].setSize(flattenedPatches.size());
        forAll(flattenedPatches, i)
            procPatches[Pstream::myProcNo()][i] = flattenedPatches[i];
        Pstream::gatherList(procPatches);
        Pstream::scatterList(procPatches);

        forAll(procPatches, procI)
        {
            if( procI == Pstream::myProcNo() )
            continue;

            const labelList& cPatches = procPatches[procI];
            label counter(0);

            while( counter < cPatches.size() )
            {
                const label patchI = cPatches[counter++];
                const label size = cPatches[counter++];
                for(label i=0;i<size;++i)
                    treatPatchesWithPatch_[patchI].appendIfNotIn
                    (
                        cPatches[counter++]
                    );
            }
        }
    }

    //- final adjusting of patches which shall be treated together
    boolList confirmed(treatPatchesWithPatch_.size(), false);
    forAll(treatPatchesWithPatch_, patchI)
    {
        if( treatPatchesWithPatch_[patchI].size() <= 1 )
        {
            confirmed[patchI] = true;
            continue;
        }

        if( confirmed[patchI] )
            continue;

        std::set<label> commonPatches;
        commonPatches.insert(patchI);

        DynList<label> front;
        front.append(patchI);
        confirmed[patchI] = true;

        while( front.size() )
        {
            const label fPatch = front.removeLastElement();

            forAll(treatPatchesWithPatch_[fPatch], i)
            {
                const label patchJ = treatPatchesWithPatch_[fPatch][i];

                if( confirmed[patchJ] )
                    continue;

                front.append(patchJ);
                confirmed[patchJ] = true;
                commonPatches.insert(patchJ);
                forAll(treatPatchesWithPatch_[patchJ], j)
                    commonPatches.insert(treatPatchesWithPatch_[patchJ][j]);
            }
        }

        forAllConstIter(std::set<label>, commonPatches, it)
        {
            const label patchJ = *it;

            treatPatchesWithPatch_[patchJ].clear();
            forAllConstIter(std::set<label>, commonPatches, iter)
                treatPatchesWithPatch_[patchJ].append(*iter);
        }
    }

    # ifdef DEBUGLayer
    for(it=edgeClassification.begin();it!=edgeClassification.end();++it)
    {
        const std::pair<label, label>& edgePair = it->first;
        const Pair<label>& nConvexAndConcave = it->second;
        Info << "Pair of patches " << edgePair.first << " "
            << edgePair.second << " is " << nConvexAndConcave << endl;
    }

    Info << "Patch names " << patchNames_ << endl;
    Info << "Treat patches with patch " << treatPatchesWithPatch_ << endl;

    label layerI(0), subsetId;
    boolList usedPatch(treatPatchesWithPatch_.size(), false);
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    forAll(treatPatchesWithPatch_, patchI)
    {
        if( usedPatch[patchI] || (boundaries[patchI].patchSize() == 0) )
            continue;

        Info << "Adding layer subset " << layerI
             << " for patch " << patchI << endl;
        usedPatch[patchI] = true;
        subsetId = mesh_.addFaceSubset("layer_"+help::scalarToText(layerI));
        ++layerI;

        forAll(treatPatchesWithPatch_[patchI], i)
        {
            const label cPatch = treatPatchesWithPatch_[patchI][i];
            usedPatch[cPatch] = true;

            label start = boundaries[cPatch].patchStart();
            const label size = boundaries[cPatch].patchSize();
            for(label i=0;i<size;++i)
                mesh_.addFaceToSubset(subsetId, start++);
        }
    }

    mesh_.write();
    # endif

    geometryAnalysed_ = true;
}

void boundaryLayers::addLayerForPatch(const label patchLabel)
{
    if( treatedPatch_[patchLabel] )
        return;

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    if( returnReduce(boundaries[patchLabel].patchSize(), sumOp<label>()) == 0 )
        return;

    boolList treatPatches(boundaries.size(), false);
    if( patchWiseLayers_ )
    {
        forAll(treatPatchesWithPatch_[patchLabel], pI)
            treatPatches[treatPatchesWithPatch_[patchLabel][pI]] = true;
    }
    else
    {
        forAll(treatedPatch_, patchI)
            if( !treatedPatch_[patchI] )
                treatPatches[patchI] = true;
    }

    newLabelForVertex_.setSize(nPoints_);
    newLabelForVertex_ = -1;
    otherVrts_.clear();
    patchKey_.clear();

    createNewVertices(treatPatches);

    createNewFacesAndCells(treatPatches);

    forAll(treatPatches, patchI)
        if( treatPatches[patchI] )
            treatedPatch_[patchI] = true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh reference
boundaryLayers::boundaryLayers
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
    msePtr_(NULL),
    meshPartitionerPtr_(NULL),
    patchWiseLayers_(true),
    terminateLayersAtConcaveEdges_(false),
    is2DMesh_(false),
    patchNames_(),
    treatedPatch_(),
    treatPatchesWithPatch_(),
    newLabelForVertex_(),
    otherVrts_(),
    patchKey_(),
    nPoints_(mesh.points().size()),
    geometryAnalysed_(false)
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    patchNames_.setSize(boundaries.size());
    forAll(boundaries, patchI)
        patchNames_[patchI] = boundaries[patchI].patchName();

    treatedPatch_.setSize(boundaries.size());
    treatedPatch_ = false;

    treatPatchesWithPatch_.setSize(boundaries.size());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

boundaryLayers::~boundaryLayers()
{
    clearOut();

    if( Pstream::parRun() )
        polyMeshGenModifier(mesh_).removeUnusedVertices();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayers::addLayerForPatch(const word& patchName)
{
    if( !geometryAnalysed_ )
        findPatchesToBeTreatedTogether();

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    forAll(boundaries, patchI)
        if( boundaries[patchI].patchName() == patchName )
            addLayerForPatch(patchI);
}

void boundaryLayers::createOTopologyLayers()
{
    patchWiseLayers_ = false;
}

void boundaryLayers::terminateLayersAtConcaveEdges()
{
    terminateLayersAtConcaveEdges_ = true;
}

void boundaryLayers::activate2DMode()
{
    polyMeshGen2DEngine mesh2DEngine(mesh_);
    const boolList& zMinPoint = mesh2DEngine.zMinPoints();
    const boolList& zMaxPoint = mesh2DEngine.zMaxPoints();

    const faceList::subList& bFaces = surfaceEngine().boundaryFaces();
    const labelList& facePatch = surfaceEngine().boundaryFacePatches();

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
            treatedPatch_[patchI] = true;
        }
    }

    forAll(treatPatchesWithPatch_, patchI)
    {
        DynList<label>& patches = treatPatchesWithPatch_[patchI];

        for(label i=patches.size()-1;i>=0;--i)
            if( treatedPatch_[patches[i]] )
                patches.removeElement(i);
    }

    is2DMesh_ = true;
}

void boundaryLayers::addLayerForAllPatches()
{
    if( !geometryAnalysed_ )
        findPatchesToBeTreatedTogether();

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    if( !patchWiseLayers_ )
    {
        forAll(boundaries, patchI)
            addLayerForPatch(patchI);
    }
    else
    {
        newLabelForVertex_.setSize(nPoints_);
        newLabelForVertex_ = -1;
        otherVrts_.clear();
        patchKey_.clear();

        //- avoid generating bnd layer at empty patches in case of 2D meshing
        label counter(0);
        forAll(treatedPatch_, patchI)
            if( !treatedPatch_[patchI] )
                ++counter;

        labelList treatedPatches(counter);
        counter = 0;
        forAll(treatedPatch_, i)
            if( !treatedPatch_[i] )
                treatedPatches[counter++] = i;

        //- create bnd layer vertices
        createNewVertices(treatedPatches);

        //- create bnd layer cells
        createLayerCells(treatedPatches);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
