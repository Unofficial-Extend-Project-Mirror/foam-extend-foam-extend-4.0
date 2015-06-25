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
#include "VRWGraph.H"
#include "Ostream.H"
#include "meshOctreeSlot.H"


// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Static data * * * * * * * * * * * * * * * * * * * //

const label meshOctreeCube::hOrder_[24][8] =
{
    {0, 1, 3, 2, 4, 5, 7, 6}, // Morton addressing Z-order
    //{6, 7, 4, 5, 1, 0, 3, 2},
    {6, 2, 1, 5, 4, 0, 3, 7},
    {6, 7, 3, 2, 1, 0, 4, 5},
    {3, 7, 6, 2, 1, 5, 4, 0},
    {4, 0, 1, 5, 6, 2, 3, 7},
    {1, 0, 4, 5, 6, 7, 3, 2},
    {3, 7, 4, 0, 1, 5, 6, 2},
    {6, 2, 3, 7, 4, 0, 1, 5},
    {3, 2, 6, 7, 4, 5, 1, 0},
    {3, 2, 1, 0, 4, 5, 6, 7},
    {6, 5, 4, 7, 3, 0, 1, 2},
    {1, 2, 6, 5, 4, 7, 3, 0},
    {3, 0, 4, 7, 6, 5, 1, 2},
    {4, 0, 3, 7, 6, 2, 1, 5},
    {1, 2, 3, 0, 4, 7, 6, 5},
    {6, 5, 1, 2, 3, 0, 4, 7},
    {1, 5, 6, 2, 3, 7, 4, 0},
    {1, 5, 4, 0, 3, 7, 6, 2},
    {4, 5, 6, 7, 3, 2, 1, 0},
    {1, 0, 3, 2, 6, 7, 4, 5},
    {3, 0, 1, 2, 6, 5, 4, 7},
    {4, 5, 1, 0, 3, 2, 6, 7},
    {4, 7, 6, 5, 1, 2, 3, 0},
    {4, 7, 3, 0, 1, 2, 6, 5}
};

const label meshOctreeCube::hOrient_[24][8] =
{
    {0, 0, 0, 0, 0, 0, 0, 0}, // Morton addressing Z-order
    //{1, 2, 0, 3, 4, 0, 5, 6},
    {0, 7, 1, 8, 5, 1, 4, 9},
    {15, 0, 2, 22, 20, 2, 19, 23},
    {20, 6, 3, 23, 15, 3, 16, 22},
    {22, 13, 4, 12, 11, 4, 1, 20},
    {11, 19, 5, 20, 22, 5, 0, 12},
    {9, 3, 6, 2, 21, 6, 17, 0},
    {10, 1, 7, 11, 12, 7, 13, 14},
    {12, 9, 8, 14, 10, 8, 18, 11},
    {6, 8, 9, 7, 17, 9, 21, 1},
    {7, 15, 10, 16, 13, 10, 12, 17},
    {5, 14, 11, 9, 0, 11, 22, 8},
    {8, 20, 12, 19, 18, 12, 10, 5},
    {18, 4, 13, 5, 8, 13, 7, 19},
    {17, 11, 14, 1, 6, 14, 23, 7},
    {2, 10, 15, 18, 19, 15, 20, 21},
    {19, 17, 16, 21, 2, 16, 3, 18},
    {14, 16, 17, 15, 23, 17, 6, 10},
    {13, 21, 18, 17, 7, 18, 8, 16},
    {16, 5, 19, 4, 3, 19, 2, 13},
    {3, 12, 20, 13, 16, 20, 15, 4},
    {23, 18, 21, 10, 14, 21, 9, 15},
    {4, 23, 22, 6, 1, 22, 11, 3},
    {21, 22, 23, 0, 9, 23, 14, 2}
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from coordinates and level
meshOctreeCube::meshOctreeCube(const meshOctreeCubeCoordinates& cc)
:
    meshOctreeCubeBasic(cc),
    activeSlotPtr_(NULL),
    subCubesPtr_(NULL),
    cubeLabel_(-1),
    containedElementsLabel_(-1),
    containedEdgesLabel_(-1)
{}

meshOctreeCube::meshOctreeCube
(
    const meshOctreeCubeCoordinates& cc,
    const label nElmts,
    meshOctreeSlot* slotPtr
)
:
    meshOctreeCubeBasic(cc),
    activeSlotPtr_(slotPtr),
    subCubesPtr_(NULL),
    cubeLabel_(0),
    containedElementsLabel_(0),
    containedEdgesLabel_(-1)
{
    slotPtr->containedTriangles_.setSize(1);
    slotPtr->containedTriangles_.setRowSize(0, nElmts),
    slotPtr->containedEdges_.setSize(0);

    for(label i=0;i<nElmts;++i)
        slotPtr->containedTriangles_(0, i) = i;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshOctreeCube::~meshOctreeCube()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FixedList<meshOctreeCube*, 8> meshOctreeCube::subCubes() const
{
    if( !subCubesPtr_ )
        FatalErrorIn
        (
            "inline  FixedList<meshOctreeCube*, 8>&"
            " meshOctreeCube::subCubes() const"
        ) << "Sub cubes do not exist!" << abort(FatalError);

    FixedList<meshOctreeCube*, 8> ret;

    for(label i=0;i<8;++i)
        ret[i] = subCube(i);

    return ret;
}

Ostream& operator<<(Ostream& os, const meshOctreeCube& oc)
{
    const meshOctreeCubeCoordinates cc(oc);
    os << cc;

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
