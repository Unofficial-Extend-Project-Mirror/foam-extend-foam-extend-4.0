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
#include "IOdictionary.H"
#include "meshOctree.H"
#include "demandDrivenData.H"

//#define DEBUGAutoRef

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeAddressing::checkGluedRegions()
{
    if( !useDATABoxes_ )
        return;

    if( meshDict_.found("checkForGluedMesh") )
    {
        if( !readBool(meshDict_.lookup("checkForGluedMesh")) )
            return;
    }
    else
    {
        return;
    }

    Info << "Removing glued regions" << endl;
    const List<direction>& boxType = this->boxType();
    const VRWGraph& nodeLabels = this->nodeLabels();
    const List<direction>& nodeType = this->nodeType();

    const VRWGraph& edgeLeaves = this->edgeLeaves();
    const VRWGraph& leafEdges = this->leafEdges();
    const FRWGraph<label, 8>& nodeLeaves = this->nodeLeaves();
    const LongList<edge>& edges = this->octreeEdges();

    DynList<label> neighbours;
    labelLongList removeBox;
    forAll(boxType, leafI)
    {
        if(
            octree_.hasContainedTriangles(leafI) &&
            (boxType[leafI] & MESHCELL)
        )
        {
            //- mark the initial INNERNODE
            labelHashSet innerNodes;
            forAllRow(nodeLabels, leafI, nodeI)
                if( nodeType[nodeLabels(leafI, nodeI)] & INNERNODE )
                {
                    innerNodes.insert(nodeLabels(leafI, nodeI));
                    break;
                }

            //- mark all INNERNODEs for which it is possible to walk along
            //- cubes edges without crossing any boundary
            bool finished;
            do
            {
                finished = true;

                forAllRow(leafEdges, leafI, leI)
                {
                    const label edgeI = leafEdges(leafI, leI);
                    const edge& edg = edges[edgeI];
                    const label s = edg[0];
                    const label e = edg[1];

                    if(
                        (nodeType[s] & INNERNODE) && (nodeType[e] & INNERNODE)
                        && (innerNodes.found(s) ^ innerNodes.found(e))
                    )
                    {
                        bool foundInside(false);
                        forAllRow(edgeLeaves, edgeI, elI)
                        {
                            const label ecLabel = edgeLeaves(edgeI, elI);
                            if( ecLabel < 0 )
                                continue;
                            if( octree_.hasContainedTriangles(ecLabel) )
                                continue;
                            if(
                                !(
                                    octree_.returnLeaf(ecLabel).cubeType() &
                                    meshOctreeCubeBasic::INSIDE
                                )
                            )
                                continue;

                            foundInside = true;
                            break;
                        }

                        if( foundInside )
                        {
                            innerNodes.insert(s);
                            innerNodes.insert(e);
                            finished = false;
                        }
                    }
                }

            } while( !finished );

            labelHashSet permissibleNeighbours;
            forAllConstIter(labelHashSet, innerNodes, it)
            {
                const label nodeI = it.key();
                forAllRow(nodeLeaves, nodeI, nlI)
                    permissibleNeighbours.insert(nodeLeaves(nodeI, nlI));
            }

            if( permissibleNeighbours.size() )
            {
                for(label i=0;i<6;++i)
                {
                    neighbours.clear();
                    octree_.findNeighboursInDirection(leafI, i, neighbours);

                    forAll(neighbours, neiI)
                    {
                        const label nei = neighbours[neiI];
                        if( nei < 0 )
                            continue;
                        if( !(boxType[nei] & MESHCELL) )
                            continue;
                        if( !permissibleNeighbours.found(nei) )
                        {
                            removeBox.append(leafI);
                            break;
                        }
                    }
                }
            }
            else
            {
                removeBox.append(leafI);
            }
        }
    }

    forAll(removeBox, i)
        (*boxTypePtr_)[removeBox[i]] = BOUNDARY;

    if( Pstream::parRun() )
    {
        LongList<meshOctreeCubeCoordinates> checkBoundary;
        forAll(removeBox, i)
        {
            const meshOctreeCubeBasic& oc = octree_.returnLeaf(removeBox[i]);
            checkBoundary.append(oc.coordinates());
        }

        LongList<meshOctreeCubeCoordinates> receivedBoundary;
        octree_.exchangeRequestsWithNeighbourProcessors
        (
            checkBoundary,
            receivedBoundary
        );

        forAll(receivedBoundary, i)
        {
            const label cLabel =
                octree_.findLeafLabelForPosition(receivedBoundary[i]);
            if( cLabel < 0 )
                continue;
            (*boxTypePtr_)[cLabel] = BOUNDARY;
        }
    }

    Info << "Finished removing glued regions" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
