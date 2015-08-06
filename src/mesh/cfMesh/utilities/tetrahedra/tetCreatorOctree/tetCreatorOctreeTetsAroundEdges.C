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

#include "tetCreatorOctree.H"
#include "demandDrivenData.H"
#include "meshOctree.H"

//#define DEBUGTets

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetCreatorOctree::createTetsAroundEdges()
{
    Info << "Creating tets around edges" << endl;

    const labelList& cubeLabel = *cubeLabelPtr_;
    const VRWGraph& nodeLabels = octreeCheck_.nodeLabels();
    const FRWGraph<label, 8>& pLeaves = octreeCheck_.nodeLeaves();
    const VRWGraph& subNodeLabels = *subNodeLabelsPtr_;
    const VRWGraph& faceCentreLabel = *faceCentreLabelPtr_;
    const meshOctree& octree = octreeCheck_.octree();
    const FixedList<FixedList<meshOctreeCubeCoordinates, 8>, 8>& vlPos =
        octree.positionsOfLeavesAtNodes();

    //- find maximum refinement level of octree leaves attached to each vertex
    List<direction> nodeLevel(octreeCheck_.numberOfNodes());

    forAll(pLeaves, nodeI)
    {
        direction level(0);

        for(label plI=0;plI<8;++plI)
        {
            const label leafI = pLeaves(nodeI, plI);

            if( leafI < 0 )
                continue;

            level = Foam::max(level, octree.returnLeaf(leafI).level());
        }

        nodeLevel[nodeI] = level;
    }

    //- start creating tets around edges which have both vertices at the same
    //- refinement level which is equal to the max refinement level of boxes
    //- incident to such edges
    forAllReverse(sortedLeaves_, levelI)
    {
        const labelLongList& curLevelLeaves = sortedLeaves_[levelI];

        forAll(curLevelLeaves, leafI)
        {
            const label curLeaf = curLevelLeaves[leafI];

            if( cubeLabel[curLeaf] == -1 )
                continue;

            const meshOctreeCubeCoordinates& cc =
                octree.returnLeaf(curLeaf).coordinates();

            # ifdef DEBUGTets
            Info << "Search cube " << curLeaf << " has coordinates "
                << octree.returnLeaf(curLeaf).coordinates() << endl;
            Info << "Node labels for cube are " << nodeLabels[curLeaf] << endl;
            # endif

            //- start checking edges
            for(label eI=0;eI<12;++eI)
            {
                const label startNode =
                    meshOctreeCubeCoordinates::edgeNodes_[eI][0];

                const label start = nodeLabels(curLeaf, startNode);
                const label end =
                    nodeLabels
                    (
                        curLeaf,
                        meshOctreeCubeCoordinates::edgeNodes_[eI][1]
                    );

                # ifdef DEBUGTets
                Info << "Creating tets around edge " << eI << endl;
                Info << "Edge nodes are " << start << " and " << end << endl;
                Info << "Coordinates start " << tetPoints_[start]
                     << " end " << tetPoints_[end] << endl;
                # endif

                bool create(true);

                if(
                    (nodeLevel[start] == direction(levelI)) &&
                    (nodeLevel[end] == direction(levelI))
                )
                {
                    //- edge has both vertices at the same refinement level
                    //- as the current leaf
                    FixedList<label, 4> edgeCubes;
                    const label fI = 2*(eI/4)+1;

                    const label* fNodes =
                        meshOctreeCubeCoordinates::faceNodes_[fI];

                    //- store octree leaves at this edge
                    //- they are all adjacent to the start point
                    for(label i=0;i<4;++i)
                        edgeCubes[i] = pLeaves(start, fNodes[i]);

                    # ifdef DEBUGTets
                    forAll(edgeCubes, i)
                    {
                        Info << "Cube " << i << " is " << edgeCubes[i] << endl;

                        if( edgeCubes[i] < 0 )
                            continue;

                        Info << "Cubes has node labels "
                             << nodeLabels[edgeCubes[i]] << endl;
                    }
                    # endif

                    DynList<label> centreNodes;

                    forAll(edgeCubes, i)
                    {
                        const label cLabel = edgeCubes[i];

                        if( (cLabel == -1) || (cubeLabel[cLabel] == -1) )
                        {
                            centreNodes.append(-1);
                            continue;
                        }

                        const meshOctreeCubeCoordinates& oc =
                            octree.returnLeaf(cLabel).coordinates();

                        # ifdef DEBUGTets
                        Info << "Edge cube " << i << " is " << oc << endl;
                        Info << "Node labels ";
                        forAllRow(nodeLabels, cLabel, k)
                            Info << nodeLabels(cLabel, k) << " ";
                        Info << endl;
                        # endif

                        if( oc.level() == direction(levelI) )
                        {
                            if( cLabel < curLeaf )
                            {
                                create = false;
                                break;
                            }

                            # ifdef DEBUGTets
                            Info << "Adding centre label " << cubeLabel[cLabel]
                                << endl;
                            # endif

                            centreNodes.append(cubeLabel[cLabel]);

                            //- adding face centre labels
                            if( faceCentreLabel.sizeOfRow(cLabel) != 0 )
                            {
                                const label helpFace = eI/4;

                                const label fcl =
                                    faceCentreLabel
                                    (
                                        cLabel,
                                        faceCentreHelper_[helpFace][i]
                                    );

                                if( fcl != -1 )
                                    centreNodes.append(fcl);
                            }

                            # ifdef DEBUGTets
                            Info << "Centre nodes after cube " << i
                                << " are " << centreNodes << endl;
                            # endif
                        }
                        else if( oc.level() < direction(levelI) )
                        {
                            # ifdef DEBUGTets
                            Info << "Edge cube " << cLabel << endl;
                            Info << "cc " << cc << endl;
                            Info << "oc " << oc << endl;
                            Info << "Adding pos "
                                 << vlPos[startNode][fNodes[i]] << endl;
                            # endif

                            const meshOctreeCubeCoordinates sc
                            (
                                cc + vlPos[startNode][fNodes[i]]
                            );

                            # ifdef DEBUGTets
                            Info << "sc " << sc << endl;
                            # endif

                            label pos(-1);

                            for(label j=0;j<8;j++)
                            {
                                if( sc == oc.refineForPosition(j) )
                                {
                                    pos = j;
                                    break;
                                }
                            }

                            if( pos == -1 )
                                FatalErrorIn
                                (
                                    "void tetCreatorOctree::"
                                    "createTetsAroundEdges()"
                                ) << "Cannot find cube position"
                                  << abort(FatalError);

                            # ifdef DEBUGTets
                            Info << "Pos " << pos << endl;
                            # endif

                            centreNodes.append(subNodeLabels(cLabel, pos));

                            # ifdef DEBUGTets
                            Info << "Centre node " << i << " is "
                                << subNodeLabels(cLabel, pos)
                                << " coordinates "
                                << tetPoints_[subNodeLabels(cLabel, pos)]
                                << endl;
                            # endif
                        }
                    }

                    //- create tets around this edge
                    if( create )
                    {
                        const label nCentres = centreNodes.size();

                        forAll(centreNodes, i)
                        {
                            if( centreNodes[i] == -1 )
                                continue;
                            if( centreNodes[(i+1)%nCentres] == -1 )
                                continue;

                            partTet tet
                            (
                                centreNodes[i],
                                centreNodes[(i+1)%nCentres],
                                start,
                                end
                            );

                            tets_.append(tet);

                            # ifdef DEBUGTets
                            Info << "Last added tet "
                                << tets_.size()-1 <<" is "
                                << tets_[tets_.size()-1] << endl;
                            # endif
                        }
                    }
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
