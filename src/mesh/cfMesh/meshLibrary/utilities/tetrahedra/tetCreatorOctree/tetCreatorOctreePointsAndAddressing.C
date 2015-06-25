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
#include "meshOctree.H"
#include "triSurface.H"

//#define DEBUGTets

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetCreatorOctree::selectElements()
{
    const List<direction>& boxType = octreeCheck_.boxType();
    const meshOctree& octree = octreeCheck_.octree();
    const boundBox& rootBox = octree.rootBox();

    //- store nodeLabels first
    const VRWGraph& nodeLabels = octreeCheck_.nodeLabels();
    tetPoints_.setSize(octreeCheck_.numberOfNodes());

    forAll(nodeLabels, leafI)
        if( nodeLabels.sizeOfRow(leafI) != 0 )
        {
            const meshOctreeCubeBasic& oc = octree.returnLeaf(leafI);
            FixedList<point, 8> lv;
            oc.vertices(rootBox, lv);

            forAll(lv, vI)
                tetPoints_[nodeLabels(leafI, vI)] = lv[vI];
        }

    //- create cubeLabel list
    if( !cubeLabelPtr_ )
        cubeLabelPtr_ = new labelList();
    labelList& cubeLabel = *cubeLabelPtr_;
    cubeLabel.setSize(octree.numberOfLeaves());
    cubeLabel = -1;

    forAll(boxType, leafI)
        if( boxType[leafI] & meshOctreeAddressing::MESHCELL )
        {
            const meshOctreeCubeBasic& oc = octree.returnLeaf(leafI);
            cubeLabel[leafI] = tetPoints_.size();
            tetPoints_.append(oc.centre(rootBox));
        }
}

void tetCreatorOctree::createPointsAndAddressing()
{
    selectElements();

    const meshOctree& octree = octreeCheck_.octree();
    const boundBox& rootBox = octree.rootBox();
    const VRWGraph& nodeLabels = octreeCheck_.nodeLabels();

    if( !subNodeLabelsPtr_ )
        subNodeLabelsPtr_ = new VRWGraph(octree.numberOfLeaves());
    VRWGraph& subNodeLabels = *subNodeLabelsPtr_;

    //- store nodeLabels
    direction maxLevel(0);
    forAll(nodeLabels, leafI)
        if( octree.returnLeaf(leafI).level() > maxLevel )
            maxLevel = octree.returnLeaf(leafI).level();

    sortedLeaves_.setSize(maxLevel+1);
    forAll(sortedLeaves_, levelI)
        sortedLeaves_[levelI].clear();

    forAll(nodeLabels, leafI)
    {
        const meshOctreeCubeBasic& oc = octree.returnLeaf(leafI);
        sortedLeaves_[oc.level()].append(leafI);
    }

    //- create subNodeLabels
    const FRWGraph<label, 8>& pointLeaves = octreeCheck_.nodeLeaves();

    forAll(pointLeaves, pointI)
    {
        bool validLeaf[8];
        direction levelI(0);
        for(label i=0;i<8;++i)
        {
            const label pointLeafI = pointLeaves(pointI, i);

            if( pointLeafI == -1 )
            {
                validLeaf[i] = false;
            }
            else
            {
                validLeaf[i] = true;
                for(label j=i+1;j<8;++j)
                    if( pointLeafI == pointLeaves(pointI, j) )
                    {
                        validLeaf[i] = false;
                        validLeaf[j] = false;
                    }

                const direction level =
                    octree.returnLeaf(pointLeafI).level();

                if( level > levelI )
                    levelI = level;
            }
        }

        for(label plI=0;plI<8;++plI)
            if( validLeaf[plI] )
            {
                const label pointLeafI = pointLeaves(pointI, plI);

                const meshOctreeCubeBasic& lc =
                    octree.returnLeaf(pointLeafI);

                if( lc.level() < levelI )
                {
                    if( subNodeLabels.sizeOfRow(pointLeafI) != 8 )
                    {
                        subNodeLabels.setRowSize(pointLeafI, 8);
                        forAllRow(subNodeLabels, pointLeafI, k)
                            subNodeLabels(pointLeafI, k) = -1;
                    }

                    subNodeLabels(pointLeafI, (7-plI)) = tetPoints_.size();
                    FixedList<point, 8> lv;
                    lc.vertices(rootBox, lv);

                    tetPoints_.append
                    (
                        0.5 * (lv[7-plI] + lc.centre(rootBox))
                    );
                }
            }
    }

    createFaceCentreLabels();

    # ifdef DEBUGTets
    forAll(nodeLabels, leafI)
    {
        forAllRow(nodeLabels, leafI, nlI)
            if( leafI != pointLeaves(nodeLabels(leafI, nlI), (7-nlI)) )
                FatalError << "Shit" << abort(FatalError);
    }
    # endif
}

void tetCreatorOctree::createFaceCentreLabels()
{
    Info << "Creating face centre labels " << endl;
    const labelList& cubeLabel = *cubeLabelPtr_;
    const VRWGraph& nodeLabels = octreeCheck_.nodeLabels();
    const FRWGraph<label, 8>& pointLeaves = octreeCheck_.nodeLeaves();
    const meshOctree& octree = octreeCheck_.octree();


    # ifdef DEBUGTets
    Info << "Node labels " << nodeLabels << endl;
    # endif

    List<direction> nodeLevel(pointLeaves.size(), direction(0));
    forAll(nodeLabels, leafI)
    {
        const direction level = octree.returnLeaf(leafI).level();

        forAllRow(nodeLabels, leafI, nlI)
        {
            const label nLabel = nodeLabels(leafI, nlI);
            # ifdef DEBUGTets
            Info << "Node label[" << leafI << "][" << nlI << "] "
                << nLabel << endl;
            # endif

            if( nodeLevel[nLabel] < level )
                nodeLevel[nLabel] = level;
        }
    }

    if( !faceCentreLabelPtr_ )
        faceCentreLabelPtr_ = new VRWGraph(cubeLabel.size());
    VRWGraph& faceCentreLabel = *faceCentreLabelPtr_;

    forAll(cubeLabel, cubeI)
        if( cubeLabel[cubeI] != -1 )
        {
            const direction level = octree.returnLeaf(cubeI).level();

            for(label i=0;i<6;++i)
            {
                if(
                    (faceCentreLabel.sizeOfRow(cubeI) != 0) &&
                    (faceCentreLabel(cubeI, i) != -1)
                )
                    continue;

                FixedList<label, 4> faceNodes;
                forAll(faceNodes, fnI)
                    faceNodes[fnI] =
                        meshOctreeCubeCoordinates::faceNodes_[i][fnI];

                label highLevelNode(-1);
                for(label j=0;j<4;++j)
                    if( nodeLevel[nodeLabels(cubeI, faceNodes[j])] > level )
                    {
                        highLevelNode = j;
                        break;
                    }

                if( highLevelNode == -1 )
                    continue;

                DynList<label> neighbours;
                octree.findNeighboursInDirection(cubeI, i, neighbours);

                if( (neighbours.size() != 1) || (neighbours[0] == -1) )
                    continue;

                if( faceCentreLabel.sizeOfRow(cubeI) == 0 )
                {
                    faceCentreLabel.setRowSize(cubeI, 6);
                    forAllRow(faceCentreLabel, cubeI, colI)
                        faceCentreLabel(cubeI, colI) = -1;
                }
                const label cNei = neighbours[0];
                if( faceCentreLabel.sizeOfRow(cNei) == 0 )
                {
                    faceCentreLabel.setRowSize(cNei, 6);
                    forAllRow(faceCentreLabel, cNei, colI)
                        faceCentreLabel(cNei, colI) = -1;
                }

                faceCentreLabel(cubeI, i) = tetPoints_.size();
                switch( i )
                {
                    case 0:
                    {
                        faceCentreLabel(cNei, 1) = tetPoints_.size();
                    } break;
                    case 1:
                    {
                        faceCentreLabel(cNei, 0) = tetPoints_.size();
                    } break;
                    case 2:
                    {
                        faceCentreLabel(cNei, 3) = tetPoints_.size();
                    } break;
                    case 3:
                    {
                        faceCentreLabel(cNei, 2) = tetPoints_.size();
                    } break;
                    case 4:
                    {
                        faceCentreLabel(cNei, 5) = tetPoints_.size();
                    } break;
                    case 5:
                    {
                        faceCentreLabel(cNei, 4) = tetPoints_.size();
                    } break;
                };

                point p(vector::zero);
                for(label j=0;j<4;++j)
                    p += tetPoints_[nodeLabels(cubeI, faceNodes[j])];
                p /= 4;
                tetPoints_.append(p);
            }
        }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
