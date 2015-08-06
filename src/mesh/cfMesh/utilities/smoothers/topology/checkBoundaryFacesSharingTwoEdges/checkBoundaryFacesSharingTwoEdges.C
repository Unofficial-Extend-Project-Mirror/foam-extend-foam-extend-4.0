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

#include "checkBoundaryFacesSharingTwoEdges.H"
#include "demandDrivenData.H"
#include "helperFunctionsPar.H"
#include "meshSurfaceEngine.H"
#include "decomposeFaces.H"
#include "decomposeCells.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGCheck

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //

void checkBoundaryFacesSharingTwoEdges::createMeshSurface() const
{
    meshSurfacePtr_ = new meshSurfaceEngine(mesh_);
}

void checkBoundaryFacesSharingTwoEdges::findFacesAtBndEdge()
{
    const meshSurfaceEngine& mse = meshSurface();

    const labelList& bp = mse.bp();
    const edgeList& edges = mse.edges();
    const VRWGraph& pointEdges = mse.boundaryPointEdges();

    const label nIntFaces = mesh_.nInternalFaces();
    const faceListPMG& faces = mesh_.faces();

    //- find the internal faces attached to the boundary points
    removeBndPoint_.setSize(pointEdges.size());
    removeBndPoint_ = true;

    # ifdef USE_OMP
    # pragma omp parallel for if( nIntFaces > 100 ) schedule(dynamic, 20)
    # endif
    for(label fI=0;fI<nIntFaces;++fI)
    {
        const face& f = faces[fI];

        forAll(f, pI)
        {
            const label bpI = bp[f[pI]];

            if( bpI < 0 )
                continue;

            if( nBndFacesAtBndPoint_[bpI] == 2 )
            {
                const edge ePrev = f.faceEdge(f.rcIndex(pI));
                const edge eNext = f.faceEdge(pI);

                bool foundNext(false), foundPrev(false);
                forAllRow(pointEdges, bpI, peI)
                {
                    const label beI = pointEdges(bpI, peI);

                    if( edges[beI] == ePrev )
                    {
                        foundPrev = true;
                    }
                    else if( edges[beI] == eNext )
                    {
                        foundNext = true;
                    }
                }

                if( !(foundPrev && foundNext) )
                {
                    removeBndPoint_[bpI] = false;
                }
            }
            else
            {
                removeBndPoint_[bpI] = false;
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- check processor faces
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label end = start + procBoundaries[patchI].patchSize();

            # ifdef USE_OMP
            # pragma omp parallel for schedule(dynamic, 10)
            # endif
            for(label faceI=start;faceI<end;++faceI)
            {
                const face& f = faces[faceI];

                forAll(f, pI)
                {
                    const label bpI = bp[f[pI]];

                    if( bpI < 0 )
                        continue;

                    if( nBndFacesAtBndPoint_[bpI] == 2 )
                    {
                        const edge ePrev = f.faceEdge(f.rcIndex(pI));
                        const edge eNext = f.faceEdge(pI);

                        bool foundNext(false), foundPrev(false);
                        forAllRow(pointEdges, bpI, peI)
                        {
                            const label beI = pointEdges(bpI, peI);

                            if( edges[beI] == ePrev )
                            {
                                foundPrev = true;
                            }
                            else if( edges[beI] == eNext )
                            {
                                foundNext = true;
                            }
                        }

                        if( !(foundPrev && foundNext) )
                            removeBndPoint_[bpI] = false;
                    }
                    else
                    {
                        removeBndPoint_[bpI] = false;
                    }
                }
            }
        }

        //- make sure that all processors have the same information
        const DynList<label>& bpNei = mse.bpNeiProcs();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

        std::map<label, labelLongList> exchangeData;
        forAll(bpNei, i)
            exchangeData.insert(std::make_pair(bpNei[i], labelLongList()));

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            if( removeBndPoint_[bpI] )
                continue;

            //- the point shall not be removed
            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);
                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append(it.key());
            }
        }

        //- exchange data
        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        //- set remove flag to false
        forAll(receivedData, i)
            removeBndPoint_[globalToLocal[receivedData[i]]] = false;
    }
}

void checkBoundaryFacesSharingTwoEdges::findBndFacesAtBndVertex()
{
    const meshSurfaceEngine& mse = meshSurface();
    const VRWGraph& pointFaces = mse.pointFaces();

    nBndFacesAtBndPoint_.setSize(pointFaces.size());
    nBndFacesAtBndPoint_ = 0;

    forAll(nBndFacesAtBndPoint_, bpI)
        nBndFacesAtBndPoint_[bpI] = pointFaces.sizeOfRow(bpI);

    if( Pstream::parRun() )
    {
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();
        const DynList<label>& neiProcs = mse.bpNeiProcs();

        //- create data that shall be exhcnaged
        std::map<label, labelLongList> exchangeData;
        forAll(neiProcs, i)
            exchangeData.insert(std::make_pair(neiProcs[i], labelLongList()));

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);
                if( neiProc == Pstream::myProcNo() )
                    continue;

                labelLongList& data = exchangeData[neiProc];
                data.append(it.key());
                data.append(nBndFacesAtBndPoint_[bpI]);
            }
        }

        //- exchange data with other processors
        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label bpI = globalToLocal[receivedData[counter++]];
            nBndFacesAtBndPoint_[bpI] += receivedData[counter++];
        }
    }
}

void checkBoundaryFacesSharingTwoEdges::removeExcessiveVertices()
{
    const labelList& bp = meshSurface().bp();
    const faceListPMG& faces = mesh_.faces();

    //- remove points which can be safely be removed
    //- internal faces
    const label nIntFaces = mesh_.nInternalFaces();

    # ifdef USE_OMP
    # pragma omp parallel for if( nIntFaces > 100 ) schedule(dynamic, 10)
    # endif
    for(label faceI=0;faceI<nIntFaces;++faceI)
    {
        const face& f = faces[faceI];

        DynList<label> newF;
        forAll(f, pI)
        {
            const label bpI = bp[f[pI]];

            if(
                (bpI >= 0) && removeBndPoint_[bpI] &&
                (nBndFacesAtBndPoint_[bpI] == 2)
            )
                continue;

            newF.append(f[pI]);
        }

        if( newF.size() < f.size() )
        {
            face& mf = const_cast<face&>(f);
            mf.setSize(newF.size());
            forAll(mf, i)
                mf[i] = newF[i];
        }
    }

    //- boundary faces
    forAll(mesh_.boundaries(), patchI)
    {
        const label start = mesh_.boundaries()[patchI].patchStart();
        const label end = start + mesh_.boundaries()[patchI].patchSize();

        # ifdef USE_OMP
        # pragma omp parallel for if( end - start > 100 ) \
        schedule(dynamic, 10)
        # endif
        for(label faceI=start;faceI<end;++faceI)
        {
            const face& f = faces[faceI];

            DynList<label> newF;
            forAll(f, pI)
            {
                const label bpI = bp[f[pI]];

                if( removeBndPoint_[bpI] && (nBndFacesAtBndPoint_[bpI] == 2) )
                    continue;

                newF.append(f[pI]);
            }

            if( newF.size() < f.size() )
            {
                face& mf = const_cast<face&>(f);
                mf.setSize(newF.size());
                forAll(mf, i)
                    mf[i] = newF[i];
            }
        }
    }

    //- processor boundaries
    forAll(mesh_.procBoundaries(), patchI)
    {
        const processorBoundaryPatch& patch = mesh_.procBoundaries()[patchI];
        const label start = patch.patchStart();
        const label end = start + patch.patchSize();

        # ifdef USE_OMP
        # pragma omp parallel for if( patch.patchSize() > 100 ) \
        schedule(dynamic, 10)
        # endif
        for(label faceI=start;faceI<end;++faceI)
        {
            const face& f = faces[faceI];

            DynList<label> newF;
            forAll(f, pI)
            {
                const label bpI = bp[f[pI]];

                if(
                    (bpI >= 0) && removeBndPoint_[bpI] &&
                    (nBndFacesAtBndPoint_[bpI] == 2)
                )
                    continue;

                newF.append(f[pI]);
            }

            if( newF.size() < f.size() )
            {
                face& mf = const_cast<face&>(f);
                mf.setSize(newF.size());

                if( !patch.owner() && (newF[0] != f[0]) )
                {
                    forAll(mf, i)
                        mf[i] = newF[mf.rcIndex(i)];
                }
                else
                {
                    forAll(mf, i)
                        mf[i] = newF[i];
                }
            }
        }
    }
}

label checkBoundaryFacesSharingTwoEdges::findBndFacesForDecomposition
(
    boolList& decomposeFace
)
{
    const meshSurfaceEngine& mse = meshSurface();
    const labelList& bp = mse.bp();
    const faceList::subList& bFaces = mse.boundaryFaces();

    label nDecomposed(0);
    const label nIntFaces = mesh_.nInternalFaces();

    # ifdef USE_OMP
    # pragma omp parallel for if( bFaces.size() > 100 ) \
    schedule(dynamic, 10) reduction(+ : nDecomposed)
    # endif
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        forAll(bf, pI)
        {
            const label bpI = bp[bf[pI]];

            if( nBndFacesAtBndPoint_[bpI] == 2 )
            {
                ++nDecomposed;
                decomposeFace[nIntFaces+bfI] = true;
            }
        }
    }

    reduce(nDecomposed, sumOp<label>());

    return nDecomposed;
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //
// Constructors

checkBoundaryFacesSharingTwoEdges::checkBoundaryFacesSharingTwoEdges
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
    meshSurfacePtr_(NULL),
    nBndFacesAtBndPoint_(),
    removeBndPoint_()
{
}

// * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * * * * //

checkBoundaryFacesSharingTwoEdges::~checkBoundaryFacesSharingTwoEdges()
{
    deleteDemandDrivenData(meshSurfacePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkBoundaryFacesSharingTwoEdges::findPoints(labelHashSet& badPoints)
{
    badPoints.clear();

    findBndFacesAtBndVertex();

    const labelList& bPoints = meshSurface().boundaryPoints();
    forAll(nBndFacesAtBndPoint_, bpI)
    {
        if( nBndFacesAtBndPoint_[bpI] != 2 )
            continue;

        badPoints.insert(bPoints[bpI]);
    }
}

bool checkBoundaryFacesSharingTwoEdges::improveTopology()
{
    bool changed(false);

    findBndFacesAtBndVertex();

    findFacesAtBndEdge();

    removeExcessiveVertices();

    boolList decomposeFace(mesh_.faces().size(), false);
    const label nDecomposed = findBndFacesForDecomposition(decomposeFace);

    Info << "Marked " << nDecomposed << " faces for decomposition" << endl;

    if( nDecomposed != 0 )
    {
        //- delete the mesh surface engine
        deleteDemandDrivenData(meshSurfacePtr_);

        //- find cells which will be decomposed
        boolList decomposeCell(mesh_.cells().size(), false);
        const labelList& owner = mesh_.owner();
        forAll(decomposeFace, faceI)
        {
            if( decomposeFace[faceI] )
                decomposeCell[owner[faceI]];
        }

        //- decompose marked faces
        decomposeFaces(mesh_).decomposeMeshFaces(decomposeFace);

        //- decompose cells
        VRWGraph pRegions(mesh_.points().size());
        decomposeCells dc(mesh_);
        dc.decomposeMesh(decomposeCell);

        changed = true;
    }

    polyMeshGenModifier(mesh_).removeUnusedVertices();

    return changed;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
