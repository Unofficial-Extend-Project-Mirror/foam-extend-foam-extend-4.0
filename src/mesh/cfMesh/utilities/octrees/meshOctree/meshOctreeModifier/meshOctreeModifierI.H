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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline const meshOctree& meshOctreeModifier::octree() const
{
    return octree_;
}

inline labelList& meshOctreeModifier::neiProcsAccess()
{
    return octree_.neiProcs_;
}

inline boundBox& meshOctreeModifier::rootBoxAccess()
{
    return octree_.rootBox_;
}

inline bool& meshOctreeModifier::isRootInitialisedAccess()
{
    return octree_.isRootInitialised_;
}

inline scalar& meshOctreeModifier::searchRangeAccess()
{
    return octree_.searchRange_;
}

inline List<Pair<meshOctreeCubeCoordinates> >&
meshOctreeModifier::neiRangeAccess()
{
    return octree_.neiRange_;
}

inline meshOctreeCube& meshOctreeModifier::initialCubeAccess()
{
    return *octree_.initialCubePtr_;
}

inline List<meshOctreeSlot>& meshOctreeModifier::dataSlotsAccess()
{
    return octree_.dataSlots_;
}

inline meshOctreeCube* meshOctreeModifier::findCubeForPosition
(
    const meshOctreeCubeCoordinates& cc
) const
{
    return octree_.findCubeForPosition(cc);
}

inline void meshOctreeModifier::findLeavesContainedInBox
(
    const boundBox& bb,
    DynList<const meshOctreeCube*, 256>& containedLeaves
) const
{
    containedLeaves.clear();
    octree_.findLeavesContainedInBox(bb, containedLeaves);
}

inline LongList<meshOctreeCube*>& meshOctreeModifier::leavesAccess()
{
    return octree_.leaves_;
}

inline void meshOctreeModifier::createListOfLeaves()
{
    octree_.leaves_.clear();

    octree_.initialCubePtr_->findLeaves(octree_.leaves_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
