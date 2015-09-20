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
#include "triSurface.H"
#include "meshOctreeSlot.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline meshOctreeCube::meshOctreeCube()
:
    meshOctreeCubeBasic(),
    activeSlotPtr_(NULL),
    subCubesPtr_(NULL),
    cubeLabel_(-1),
    containedElementsLabel_(-1),
    containedEdgesLabel_(-1)
{}

inline meshOctreeCube::meshOctreeCube(const meshOctreeCube& moc)
:
    meshOctreeCubeBasic(moc.coordinates(), moc.cubeType(), moc.procNo()),
    activeSlotPtr_(moc.activeSlotPtr_),
    subCubesPtr_(moc.subCubesPtr_),
    cubeLabel_(moc.cubeLabel_),
    containedElementsLabel_(moc.containedElementsLabel_),
    containedEdgesLabel_(moc.containedEdgesLabel_)
{}

inline const meshOctreeSlot* meshOctreeCube::slotPtr() const
{
    return activeSlotPtr_;
}

inline bool meshOctreeCube::isLeaf() const
{
    if( !subCubesPtr_ )
        return true;

    return false;
}

inline meshOctreeCube* meshOctreeCube::subCube(const label scI) const
{
    return subCubesPtr_[scI];
}

inline label meshOctreeCube::cubeLabel() const
{
    return cubeLabel_;
}

inline bool meshOctreeCube::hasContainedElements() const
{
    if( containedElementsLabel_ == -1 )
        return false;

    return true;
}

inline label meshOctreeCube::containedElements() const
{
    return containedElementsLabel_;
}

inline bool meshOctreeCube::hasContainedEdges() const
{
    if( containedEdgesLabel_ == -1 )
        return false;;

    return true;
}

inline label meshOctreeCube::containedEdges() const
{
    return containedEdgesLabel_;
}

inline void meshOctreeCube::operator=(const meshOctreeCube& moc)
{
    static_cast<meshOctreeCubeCoordinates&>(*this) = moc.coordinates();
    activeSlotPtr_ = moc.activeSlotPtr_;
    subCubesPtr_ = moc.subCubesPtr_;
    cubeLabel_ = moc.cubeLabel_;
    containedElementsLabel_ = moc.containedElementsLabel_;
    containedEdgesLabel_ = moc.containedEdgesLabel_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
