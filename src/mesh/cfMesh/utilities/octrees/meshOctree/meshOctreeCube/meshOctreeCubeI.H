/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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
    activeSlotPtr_(nullptr),
    subCubesPtr_(nullptr),
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
