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

#include "meshOctreeModifier.H"
#include "helperFunctionsPar.H"
#include "triSurf.H"

#include <map>

//#define OCTREE_DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeModifier::refineTreeForCoordinates
(
    const meshOctreeCubeCoordinates& cc,
    const short procNo,
    const direction cubeType
)
{
    const label cpx = cc.posX();
    const label cpy = cc.posY();
    const label cpz = cc.posZ();
    const direction l = cc.level();

    # ifdef OCTREE_DEBUG
    const label levelLimiter = (1 << l);
    if(
        (cpx >= levelLimiter) || (cpx < 0) ||
        (cpy >= levelLimiter) || (cpy < 0) ||
        (cpz >= levelLimiter) || (cpz < 0)
    )
    {
        FatalErrorIn
        (
            "void meshOctree::refineTreeForCoordinates("
            "const meshOctreeCubeCoordinates& cc)"
        ) << "Trying to add an invalid cube!" << abort(FatalError);
    }
    # endif

    meshOctreeCube* nei(octree_.initialCubePtr_);

    for(label i=(l-1);i>=0;--i)
    {
        const label levelLimiter = (1 << i);

        label scI(0);

        if( cpx & levelLimiter )
            scI |= 1;
        if( cpy & levelLimiter )
            scI |= 2;
        if( cpz & levelLimiter )
            scI |= 4;

        if( nei->isLeaf() )
        {
            //- refine the missing cube
            nei->refineMissingCube
            (
                octree_.surface_,
                octree_.rootBox_,
                scI
            );

            nei = nei->subCube(scI);
        }
        else
        {
            if( !nei->subCube(scI) )
            {
                //- create the needed cube if it is not present
                nei->refineMissingCube
                (
                    octree_.surface_,
                    octree_.rootBox_,
                    scI
                );
            }

            nei = nei->subCube(scI);
        }
    }

    nei->setProcNo(procNo);
    nei->setCubeType(cubeType);
}

void meshOctreeModifier::refineTreeForCoordinates
(
    const meshOctreeCubeCoordinates& cc,
    const labelList& containedTriangles,
    const labelList& containedEdges,
    const short procNo,
    const direction cubeType
)
{
    const label cpx = cc.posX();
    const label cpy = cc.posY();
    const label cpz = cc.posZ();
    const direction l = cc.level();

    # ifdef OCTREE_DEBUG
    const label levelLimiter = (1 << l);
    if(
        (cpx >= levelLimiter) || (cpx < 0) ||
        (cpy >= levelLimiter) || (cpy < 0) ||
        (cpz >= levelLimiter) || (cpz < 0)
    )
    {
        FatalErrorIn
        (
            "void meshOctree::refineTreeForCoordinates("
            "const meshOctreeCubeCoordinates& cc)"
        ) << "Trying to add an invalid cube!" << abort(FatalError);
    }
    # endif

    meshOctreeCube* nei(octree_.initialCubePtr_);

    for(label i=(l-1);i>=0;--i)
    {
        const label levelLimiter = (1 << i);

        label scI(0);

        if( cpx & levelLimiter )
            scI |= 1;
        if( cpy & levelLimiter )
            scI |= 2;
        if( cpz & levelLimiter )
            scI |= 4;

        if( nei->isLeaf() )
        {
            //- refine the missing cube
            //nei->refineMissingCube(scI, containedTrianglesI, containedEdgesI);
            nei->refineMissingCube
            (
                octree_.surface_,
                octree_.rootBox_,
                scI
            );

            nei = nei->subCube(scI);
        }
        else
        {
            meshOctreeCube* scPtr = nei->subCube(scI);

            if( !scPtr )
            {
                //- create the needed cube if it is not present
                nei->refineMissingCube
                (
                    octree_.surface_,
                    octree_.rootBox_,
                    scI
                );
            }

            nei = scPtr;
        }
    }

    nei->setProcNo(procNo);
    nei->setCubeType(cubeType);
}

void meshOctreeModifier::addLayerFromNeighbouringProcessors()
{
    if( !Pstream::parRun() )
        return;

    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;
    const labelList& neiProcs = octree_.neiProcs_;
    const List<Pair<meshOctreeCubeCoordinates> >& neiRange = octree_.neiRange_;

    forAll(leaves, leafI)
        if( leaves[leafI]->procNo() != Pstream::myProcNo() )
            return;

    Info << "Adding an additional layer of cells" << endl;

    meshOctreeCubeCoordinates minCoord, maxCoord;
    std::map<label, LongList<meshOctreeCubeBasic> > toProcs;
    forAll(neiProcs, i)
        toProcs.insert
        (
            std::make_pair(neiProcs[i], LongList<meshOctreeCubeBasic>())
        );

    //- fill the data into into the map
    forAll(leaves, leafI)
    {
        leaves[leafI]->neighbourRange(minCoord, maxCoord);

        forAll(neiProcs, procI)
        {
            if(
                (maxCoord >= neiRange[procI].first()) &&
                (minCoord <= neiRange[procI].second())
            )
                toProcs[neiProcs[procI]].append(*leaves[leafI]);
        }
    }

    # ifdef OCTREE_DEBUG
    for(label i=0;i<Pstream::nProcs();++i)
    {
        if( i == Pstream::myProcNo() )
        {
            std::map<label, LongList<meshOctreeCubeBasic> >::iterator it;
            for(it=toProcs.begin();it!=toProcs.end();++it)
            {
                Pout << "Sending " << it->second.size() << " cubes to proc "
                    << it->first << endl;
            }
        }

        returnReduce(i, sumOp<label>());
    }
    # endif

    //- exchange data with other processors
    LongList<meshOctreeCubeBasic> receivedCoordinates;
    help::exchangeMap(toProcs, receivedCoordinates, Pstream::blocking);

    # ifdef OCTREE_DEBUG
    Pout << "Received " << receivedCoordinates.size()
        << " from other procs" << endl;
    # endif

    //- cubes which share a common, face, edge or vertex are added into
    //- the current processor's tree
    DynList<label> neighbours;
    forAll(receivedCoordinates, ccI)
    {
        octree_.findAllLeafNeighbours
        (
            receivedCoordinates[ccI].coordinates(),
            neighbours
        );

        forAll(neighbours, neiI)
        {
            const label nei = neighbours[neiI];

            if( nei < 0 )
                continue;

            if( leaves[nei]->procNo() == Pstream::myProcNo() )
            {
                refineTreeForCoordinates
                (
                    receivedCoordinates[ccI].coordinates(),
                    receivedCoordinates[ccI].procNo(),
                    receivedCoordinates[ccI].cubeType()
                );

                break;
            }
        }
    }

    //- recalculate leaves
    createListOfLeaves();

    # ifdef OCTREE_DEBUG
    for(label leafI=0;leafI<octree_.numberOfLeaves();++leafI)
    {
        const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);

        if( oc.procNo() == Pstream::myProcNo() )
            continue;

        DynList<label> neighbours;
        octree_.findAllLeafNeighbours(leafI, neighbours);

        bool found(false);
        forAll(neighbours, i)
        {
            const label neiLeaf = neighbours[i];

            if( neiLeaf < 0 )
                continue;
            if( octree_.returnLeaf(neiLeaf).procNo() == Pstream::myProcNo() )
            {
                found = true;
                break;
            }
        }

        if( !found )
            FatalError << "Leaf " << leafI << " with coordinates "
                << octree_.returnLeaf(leafI)
                << " has no neighbour local at this processor"
                << abort(FatalError);
    }
    # endif

    Info << "Finished adding an additional layer of octree cubes" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
