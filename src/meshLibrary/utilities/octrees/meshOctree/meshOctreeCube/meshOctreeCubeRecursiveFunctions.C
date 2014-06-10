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

#include "meshOctreeCube.H"
#include "demandDrivenData.H"
#include "Ostream.H"
#include "meshOctree.H"

//#define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeCube::leavesInBox
(
    const boundBox& rootBox,
    const boundBox& searchingBox,
    DynList<const meshOctreeCube*, 256>& leaves
) const
{
    point min, max;
    this->cubeBox(rootBox, min, max);
    const boundBox cubeBox(min, max);

    if( cubeBox.overlaps(searchingBox) )
    {
        if( this->isLeaf() )
        {
            leaves.append(this);
        }
        else
        {
            for(label scI=0;scI<8;++scI)
            {
                meshOctreeCube* scPtr = subCubesPtr_[scI];

                if( scPtr )
                {
                    scPtr->leavesInBox
                    (
                        rootBox,
                        searchingBox,
                        leaves
                    );
                }
                else if( Pstream::parRun() )
                {
                    meshOctreeCubeCoordinates cc = refineForPosition(scI);
                    cc.cubeBox(rootBox, min, max);
                    const boundBox bb(min, max);
                    if( bb.overlaps(searchingBox) )
                        leaves.append(this);
                }
            }
        }
    }
}

void meshOctreeCube::findLeaves(LongList<meshOctreeCube*>& leaves) const
{
    if( this->isLeaf() )
    {
        meshOctreeCube* oc = const_cast<meshOctreeCube*>(this);
        cubeLabel_ = leaves.size();
        leaves.append(oc);
    }
    else
    {
        cubeLabel_ = -1;

        for(label scI=0;scI<8;++scI)
        {
            const meshOctreeCube* scPtr = subCubesPtr_[scI];

            if( scPtr )
                scPtr->findLeaves(leaves);
        }
    }
}

void meshOctreeCube::findCoordinatesOfMissingCubes
(
    LongList<meshOctreeCubeCoordinates>& coordinates
) const
{
    if( this->isLeaf() )
        return;

    for(label scI=0;scI<8;++scI)
    {
        meshOctreeCube* scPtr = subCubesPtr_[scI];

        if( scPtr )
        {
            scPtr->findCoordinatesOfMissingCubes(coordinates);
        }
        else
        {
            coordinates.append(this->refineForPosition(scI));
        }
    }
}

void meshOctreeCube::countChildCubes(label& counter) const
{
    ++counter;

    if( !this->isLeaf() )
    {
        for(label scI=0;scI<8;++scI)
        {
            meshOctreeCube* scPtr = subCubesPtr_[scI];

            if( scPtr )
            {
                scPtr->countChildCubes(counter);
            }
        }
    }
}

bool meshOctreeCube::purgeProcessorCubes(const short procNo)
{
    if( this->procNo() == ALLPROCS )
    {
        this->setProcNo(procNo);
    }

    if( this->isLeaf() )
    {
        if( this->procNo() != procNo )
        {
            return true;
        }

        return false;
    }
    else
    {
        label mergedSubcubes = 0;
        for(label scI=0;scI<8;++scI)
        {
            if( !subCubesPtr_[scI] )
            {
                mergedSubcubes |= 1 << scI;
                continue;
            }

            if( subCubesPtr_[scI]->purgeProcessorCubes(procNo) )
            {
                subCubesPtr_[scI] = NULL;
                mergedSubcubes |= 1 << scI;
            }
        }

        if( mergedSubcubes == 255 )
        {
            subCubesPtr_ = NULL;

            return true;
        }
        else
        {
            return false;
        }
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
