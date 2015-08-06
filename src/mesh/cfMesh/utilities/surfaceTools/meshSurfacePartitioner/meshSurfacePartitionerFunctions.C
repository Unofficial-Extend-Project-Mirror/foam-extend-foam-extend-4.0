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

#include "meshSurfacePartitioner.H"
#include "helperFunctionsPar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfacePartitioner::calculateCornersEdgesAndAddressing()
{
    const labelList& bPoints = meshSurface_.boundaryPoints();
    const labelList& bp = meshSurface_.bp();
    const edgeList& edges = meshSurface_.edges();
    const VRWGraph& edgeFaces = meshSurface_.edgeFaces();
    const VRWGraph& pointFaces = meshSurface_.pointFaces();

    corners_.clear();
    edgePoints_.clear();

    //- count the number of patches
    label nPatches(0);
    # ifdef USE_OMP
    # pragma omp parallel
    {
        label localMax(0);

        forAll(facePatch_, bfI)
            localMax = Foam::max(localMax, facePatch_[bfI]);

        # pragma omp critical
        nPatches = Foam::max(localMax, nPatches);
    }
    # else
    forAll(facePatch_, bfI)
        nPatches = Foam::max(nPatches, facePatch_[bfI]);
    # endif
    ++nPatches;

    //- set the size and starting creating addressing
    patchPatches_.setSize(nPatches);

    nEdgesAtPoint_.clear();
    nEdgesAtPoint_.setSize(bPoints.size());
    nEdgesAtPoint_ = 0;

    featureEdges_.clear();

    forAll(edgeFaces, edgeI)
    {
        if( edgeFaces.sizeOfRow(edgeI) != 2 )
            continue;

        const label patch0 = facePatch_[edgeFaces(edgeI, 0)];
        const label patch1 = facePatch_[edgeFaces(edgeI, 1)];

        if( patch0 != patch1 )
        {
            const edge& e = edges[edgeI];
            ++nEdgesAtPoint_[bp[e.start()]];
            ++nEdgesAtPoint_[bp[e.end()]];

            patchPatches_[patch0].insert(patch1);
            patchPatches_[patch1].insert(patch0);

            featureEdges_.insert(edgeI);
        }
    }

    if( Pstream::parRun() )
    {
        const Map<label>& otherFaceAtProc = meshSurface_.otherEdgeFaceAtProc();

        //- find patches on other procs sharing surface edges
        Map<label> otherFacePatch;

        const DynList<label>& beNeiProcs = meshSurface_.beNeiProcs();
        const Map<label>& globalToLocalEdges =
            meshSurface_.globalToLocalBndEdgeAddressing();

        std::map<label, labelLongList> exchangeData;
        forAll(beNeiProcs, i)
            exchangeData.insert(std::make_pair(beNeiProcs[i], labelLongList()));

        forAllConstIter(Map<label>, globalToLocalEdges, it)
        {
            const label beI = it();

            if( edgeFaces.sizeOfRow(beI) == 1 )
            {
                labelLongList& data = exchangeData[otherFaceAtProc[beI]];

                data.append(it.key());
                data.append(facePatch_[edgeFaces(beI, 0)]);
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        for(label i=0;i<receivedData.size();)
        {
            const label geI = receivedData[i++];
            const label patchI = receivedData[i++];

            otherFacePatch.insert(globalToLocalEdges[geI], patchI);
        }

        //- take into account feature edges at processor boundaries
        forAllConstIter(Map<label>, otherFaceAtProc, it)
        {
            const label beI = it.key();

            if( it() <= Pstream::myProcNo() )
                continue;
            if( otherFacePatch[beI] != facePatch_[edgeFaces(beI, 0)] )
            {
                const edge& e = edges[beI];
                ++nEdgesAtPoint_[bp[e.start()]];
                ++nEdgesAtPoint_[bp[e.end()]];
            }
        }

        //- gather data on all processors
        exchangeData.clear();
        const DynList<label>& bpNeiProcs = meshSurface_.bpNeiProcs();
        forAll(bpNeiProcs, i)
            exchangeData.insert(std::make_pair(bpNeiProcs[i], labelLongList()));

        const Map<label>& globalToLocal =
            meshSurface_.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = meshSurface_.bpAtProcs();
        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            forAllRow(bpAtProcs, bpI, i)
            {
                const label procI = bpAtProcs(bpI, i);

                if( procI == Pstream::myProcNo() )
                    continue;

                labelLongList& dts = exchangeData[procI];

                //- exchange data as follows:
                //- 1. global point label
                //- 2. number of feature edges connected to the vertex
                dts.append(it.key());
                dts.append(nEdgesAtPoint_[bpI]);
            }
        }

        //- exchange information
        receivedData.clear();
        help::exchangeMap(exchangeData, receivedData);

        //- add the edges from other processors to the points
        label counter(0);
        while( counter < receivedData.size() )
        {
            const label bpI = globalToLocal[receivedData[counter++]];
            const label nEdges = receivedData[counter++];

            nEdgesAtPoint_[bpI] += nEdges;
        }
    }

    //- mark edges and corners
    forAll(nEdgesAtPoint_, bpI)
    {
        if( nEdgesAtPoint_[bpI] > 2 )
        {
            corners_.insert(bpI);
        }
        else if( nEdgesAtPoint_[bpI] == 2 )
        {
            edgePoints_.insert(bpI);
        }
    }

    //- find patches at a surface points
    pointPatches_.setSize(pointFaces.size());
    forAll(pointFaces, bpI)
    {
        forAllRow(pointFaces, bpI, pfI)
            pointPatches_.appendIfNotIn(bpI, facePatch_[pointFaces(bpI, pfI)]);
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal =
            meshSurface_.globalToLocalBndPointAddressing();
        const DynList<label>& bpNeiProcs = meshSurface_.bpNeiProcs();
        const VRWGraph& bpAtProcs = meshSurface_.bpAtProcs();

        std::map<label, labelLongList> exchangeData;
        forAll(bpNeiProcs, i)
            exchangeData[bpNeiProcs[i]].clear();

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
                data.append(pointPatches_.sizeOfRow(bpI));
                forAllRow(pointPatches_, bpI, i)
                    data.append(pointPatches_(bpI, i));
            }
        }

        //- exchange data with other prcessors
        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label bpI = globalToLocal[receivedData[counter++]];
            const label size = receivedData[counter++];

            for(label i=0;i<size;++i)
                pointPatches_.appendIfNotIn(bpI, receivedData[counter++]);
        }
    }

    label counter = corners_.size();
    reduce(counter, sumOp<label>());
    Info << "Found " << counter
        << " corners at the surface of the volume mesh" << endl;
    counter = edgePoints_.size();
    reduce(counter, sumOp<label>());
    Info << "Found " << counter
        << " edge points at the surface of the volume mesh" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
