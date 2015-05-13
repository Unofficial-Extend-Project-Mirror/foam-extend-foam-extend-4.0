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

#include "meshOctree.H"
#include "boundBox.H"

//#define OCTREE_DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private member functions

label meshOctree::findLeafContainingVertex
(
    const point& p
) const
{
    # ifdef OCTREE_DEBUG
    Info << "Finding leaf for vertex " << p << endl;
    # endif

    const meshOctreeCube* ocPtr = initialCubePtr_;

    if( !ocPtr->isVertexInside(rootBox_, p) )
    {
        # ifdef OCTREE_DEBUG
        Info << "Vertex " << p << " is not in the initial cube" << endl;
        # endif
        return -1;
    }

    bool finished(false);

    do
    {
        if( ocPtr && !ocPtr->isLeaf() )
        {
            //- find a subCube containing the vertex;
            const point c = ocPtr->centre(rootBox_);

            label scI(0);

            if( p.x() >= c.x() )
                scI |= 1;
            if( p.y() >= c.y() )
                scI |= 2;
            if( !isQuadtree_ && p.z() >= c.z() )
                scI |= 4;

            ocPtr = ocPtr->subCube(scI);
        }
        else
        {
            finished = true;
        }
    } while( !finished );

    if( ocPtr )
    {
        return ocPtr->cubeLabel();
    }

    return meshOctreeCubeBasic::OTHERPROC;
}

label meshOctree::findNeighbourOverNode
(
    const meshOctreeCubeCoordinates& cc,
    const label nodeI
) const
{
    if( isQuadtree_ )
        return -1;

    const meshOctreeCubeCoordinates nc(cc + regularityPositions_[18+nodeI]);

    const meshOctreeCube* neiPtr = findCubeForPosition(nc);

    if( !neiPtr )
    {
        const label levelLimiter = (1 << cc.level());
        if(
            (nc.posX() >= levelLimiter) || (nc.posX() < 0) ||
            (nc.posY() >= levelLimiter) || (nc.posY() < 0) ||
            (!isQuadtree_ && (nc.posZ() >= levelLimiter || nc.posZ() < 0)) ||
            (isQuadtree_ && (nc.posZ() != initialCubePtr_->posZ()))
        )
        {
            return -1;
        }
        else if( Pstream::parRun() )
        {
            return meshOctreeCubeBasic::OTHERPROC;
        }

        return -1;
    }
    else if( neiPtr->isLeaf() )
    {
        # ifdef OCTREE_DEBUG
        if( leaves_[neiPtr->cubeLabel()] != neiPtr )
            FatalError << "Cube does not correspond to itself"
                << abort(FatalError);
        # endif
        return neiPtr->cubeLabel();
    }
    else
    {
        FixedList<label, 8> sc(-1);
        for(label scI=0;scI<8;++scI)
        {
            meshOctreeCube* scPtr = neiPtr->subCube(scI);
            if( scPtr )
            {
                sc[scI] = scPtr->cubeLabel();
            }
            else if( Pstream::parRun() )
            {
                sc[scI] = meshOctreeCubeBasic::OTHERPROC;
            }
        }

        return sc[7-nodeI];
    }

    FatalErrorIn
    (
        "label meshOctree::findNeighbourOverNode("
        "const meshOctreeCubeCoordinates& cc,"
        "const label nodeI) const"
    ) << "Should not be here!" << abort(FatalError);

    return -1;
}

void meshOctree::findNeighboursOverEdge
(
    const meshOctreeCubeCoordinates& cc,
    const label eI,
    DynList<label>& neighbourLeaves
) const
{
    if( isQuadtree_ && (eI >= 8) )
    {
        neighbourLeaves.append(-1);
        return;
    }

    const meshOctreeCubeCoordinates nc(cc + regularityPositions_[6+eI]);

    const meshOctreeCube* neiPtr = findCubeForPosition(nc);

    if( !neiPtr )
    {
        const label levelLimiter = (1 << cc.level());
        if(
            (nc.posX() >= levelLimiter) || (nc.posX() < 0) ||
            (nc.posY() >= levelLimiter) || (nc.posY() < 0) ||
            (!isQuadtree_ && (nc.posZ() >= levelLimiter || nc.posZ() < 0)) ||
            (isQuadtree_ && (nc.posZ() != initialCubePtr_->posZ()))
        )
        {
            neighbourLeaves.append(-1);
        }
        else if( Pstream::parRun() )
        {
            neighbourLeaves.append(meshOctreeCubeBasic::OTHERPROC);
        }
    }
    else if( neiPtr->isLeaf() )
    {
        # ifdef OCTREE_DEBUG
        if( leaves_[neiPtr->cubeLabel()] != neiPtr )
            FatalError << "Cube does not correspond to itself"
                << abort(FatalError);
        # endif
        neighbourLeaves.append(neiPtr->cubeLabel());
    }
    else
    {
        FixedList<label, 8> sc(-1);
        for(label scI=0;scI<8;++scI)
        {
            meshOctreeCube* scPtr = neiPtr->subCube(scI);

            if( scPtr )
            {
                sc[scI] = scPtr->cubeLabel();
            }
            else if( Pstream::parRun() )
            {
                sc[scI] = meshOctreeCubeBasic::OTHERPROC;
            }
        }

        const label* eNodes = meshOctreeCubeCoordinates::edgeNodes_[eI];

        if( !isQuadtree_)
        {
            neighbourLeaves.append(sc[7-eNodes[1]]);
            neighbourLeaves.append(sc[7-eNodes[0]]);
        }
        else
        {
            if( sc[7-eNodes[1]] >= 0 )
                neighbourLeaves.append(sc[7-eNodes[1]]);
            if( (sc[7-eNodes[0]] >= 0) && (sc[7-eNodes[0]] != sc[7-eNodes[1]]) )
                neighbourLeaves.append(sc[7-eNodes[0]]);
        }
    }
}

void meshOctree::findNeighboursInDirection
(
    const meshOctreeCubeCoordinates& cc,
    const label dir,
    DynList<label>& neighbourLeaves
) const
{
    if( isQuadtree_ && dir >= 4 )
    {
        neighbourLeaves.append(-1);
        return;
    }

    label cpx = cc.posX();
    label cpy = cc.posY();
    label cpz = cc.posZ();
    switch( dir )
    {
        case 0:
        {
            cpx -= 1;
        }
        break;
        case 1:
        {
            cpx += 1;
        }
        break;
        case 2:
        {
            cpy -= 1;
        }
        break;
        case 3:
        {
            cpy += 1;
        }
        break;
        case 4:
        {
            cpz -= 1;
        }
        break;
        case 5:
        {
            cpz += 1;
        }
        break;
    }

    const meshOctreeCube* neiPtr =
        findCubeForPosition
        (
            meshOctreeCubeCoordinates(cpx, cpy, cpz, cc.level())
        );

    if( !neiPtr )
    {
        const label levelLimiter = (1 << cc.level());
        if(
            (cpx >= levelLimiter) || (cpx < 0) ||
            (cpy >= levelLimiter) || (cpy < 0) ||
            (!isQuadtree_ && (cpz >= levelLimiter || cpz < 0)) ||
            (isQuadtree_ && (cpz != initialCubePtr_->posZ()))
        )
        {
            neighbourLeaves.append(-1);
        }
        else if( Pstream::parRun() )
        {
            neighbourLeaves.append(meshOctreeCubeBasic::OTHERPROC);
        }
    }
    else if( neiPtr->isLeaf() )
    {
        # ifdef OCTREE_DEBUG
        if( leaves_[neiPtr->cubeLabel()] != neiPtr )
            FatalError << "Cube does not correspond to itself"
                << abort(FatalError);
        # endif
        neighbourLeaves.append(neiPtr->cubeLabel());
    }
    else
    {
        FixedList<label, 8> sc(-1);
        for(label scI=0;scI<8;++scI)
        {
            meshOctreeCube* scPtr = neiPtr->subCube(scI);

            if( scPtr )
            {
                sc[scI] = scPtr->cubeLabel();
            }
            else if( Pstream::parRun() )
            {
                sc[scI] = meshOctreeCubeBasic::OTHERPROC;
            }
        }

        const label* fNodes = meshOctreeCubeCoordinates::faceNodes_[dir];
        for(label i=0;i<4;++i)
        {
            if( isQuadtree_ && sc[7-fNodes[i]] < 0 )
                continue;

            neighbourLeaves.append(sc[7-fNodes[i]]);
        }
    }
}

void meshOctree::findNeighboursForLeaf
(
    const meshOctreeCubeCoordinates& cc,
    DynList<label>& neighbourLeaves
) const
{
    neighbourLeaves.clear();

    const label nCubeFaces = isQuadtree_?4:6;
    for(label i=0;i<nCubeFaces;++i)
    {
        findNeighboursInDirection(cc, i, neighbourLeaves);
    }
}

void meshOctree::findAllLeafNeighbours
(
    const meshOctreeCubeCoordinates& cc,
    DynList<label>& neighbourLeaves
) const
{
    neighbourLeaves.clear();

    //- neighbours over nodes
    if( !isQuadtree_ )
    {
        for(label i=0;i<8;++i)
            neighbourLeaves.append(findNeighbourOverNode(cc, i));
    }

    //- neighbours over edges
    const label nCubeEdges = isQuadtree_?8:12;
    for(label i=0;i<nCubeEdges;++i)
        findNeighboursOverEdge(cc, i, neighbourLeaves);

    //- neighbours over faces
    const label nCubeFaces = isQuadtree_?4:6;
    for(label i=0;i<nCubeFaces;++i)
        findNeighboursInDirection(cc, i, neighbourLeaves);
}

void meshOctree::findLeavesForCubeVertex
(
    const label leafI,
    const direction vrtI,
    FixedList<label, 8>& neighbours
) const
{
    const meshOctreeCube* oc = leaves_[leafI];
    const meshOctreeCubeCoordinates cc = oc->refineForPosition(vrtI);

    FixedList<meshOctreeCubeCoordinates, 8> positions;

    for(label i=0;i<8;++i)
    {
        positions[i] = cc + vrtLeavesPos_[vrtI][i];
    }

    forAll(positions, posI)
    {
        neighbours[posI] = -1;

        const label nei = findLeafLabelForPosition(positions[posI]);

        if( (nei > -1) && leaves_[nei]->isLeaf() )
            neighbours[posI] = nei;
    }

    # ifdef OCTREE_DEBUG
    Info << "Cube " << *oc << endl;
    Info << "Vertex " << short(vrtI) << endl;
    Info << "Neighbours " << endl;
    forAll(neighbours, nI)
        if( neighbours[nI] )
            Info << *neighbours[nI] << endl;
    # endif
}

meshOctreeCube* meshOctree::findCubeForPosition
(
    const meshOctreeCubeCoordinates& cc
) const
{
    # ifdef OCTREE_DEBUG
    Info << "Finding position " << cc << endl;
    # endif

    const label cpx = cc.posX();
    const label cpy = cc.posY();
    const label cpz = cc.posZ();
    const direction l = cc.level();

    label levelLimiter = (1 << l);
    if(
        (cpx >= levelLimiter) || (cpx < 0) ||
        (cpy >= levelLimiter) || (cpy < 0) ||
        (!isQuadtree_ && (cpz >= levelLimiter || cpz < 0)) ||
        (isQuadtree_ && (cpz != initialCubePtr_->posZ()))
    )
    {
        return NULL;
    }

    meshOctreeCube* neiPtr(initialCubePtr_);

    for(label i=(l-1);i>=0;--i)
    {
        if( neiPtr && !neiPtr->isLeaf() )
        {
            levelLimiter = (1 << i);

            label scI(0);

            if( cpx & levelLimiter )
                scI |= 1;
            if( cpy & levelLimiter )
                scI |= 2;
            if( !isQuadtree_ && (cpz & levelLimiter) )
                scI |= 4;

            neiPtr = neiPtr->subCube(scI);
        }
        else
        {
            break;
        }
    }

    # ifdef OCTREE_DEBUG
    Info << "Found position is " << *neiPtr << endl;
    # endif

    return neiPtr;
}

label meshOctree::findLeafLabelForPosition
(
    const meshOctreeCubeCoordinates& cc
) const
{
    const meshOctreeCube* ocPtr = findCubeForPosition(cc);

    if( ocPtr && ocPtr->isLeaf() )
    {
        return ocPtr->cubeLabel();
    }
    else if( !ocPtr && (neiProcs_.size() != 0) )
    {
        const label levelLimiter = (1 << cc.level());
        if(
            (cc.posX() < levelLimiter) && (cc.posX() >= 0) &&
            (cc.posY() < levelLimiter) && (cc.posY() >= 0) &&
            ((!isQuadtree_ && (cc.posZ() < levelLimiter && cc.posZ() >= 0)) ||
            (isQuadtree_ && (cc.posZ() == initialCubePtr_->posZ())))
        )
        {
            return meshOctreeCubeBasic::OTHERPROC;
        }
    }

    return -1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
