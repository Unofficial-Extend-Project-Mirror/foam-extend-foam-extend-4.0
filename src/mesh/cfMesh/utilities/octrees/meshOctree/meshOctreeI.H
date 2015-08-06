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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline bool meshOctree::isQuadtree() const
{
    return isQuadtree_;
}

//- return octant vectors
inline const FixedList<Vector<label>, 8>& meshOctree::octantVectors() const
{
    return octantVectors_;
}

//- return leaves of the octree
inline label meshOctree::numberOfLeaves() const
{
    if( leaves_.size() == 0 )
        FatalErrorIn
        (
            "inline const label& meshOctree::numberOfLeaves() const"
        ) << "Leaves are not yet created!!" << exit(FatalError);

    return leaves_.size();
}

inline const meshOctreeCubeBasic& meshOctree::returnLeaf
(
    const label leafI
) const
{
    return *leaves_[leafI];
}

inline short meshOctree::leafAtProc(const label leafI) const
{
    return leaves_[leafI]->procNo();
}

inline bool meshOctree::hasContainedTriangles(const label leafI) const
{
    if( leaves_[leafI]->hasContainedElements() )
        return true;

    return false;
}

inline void meshOctree::containedTriangles
(
    const label leafI,
    DynList<label>& triangles
) const
{
    triangles.clear();

    if( !leaves_[leafI]->hasContainedElements() )
    {
        triangles.clear();
        return;
    }

    const VRWGraph& containedTriangles =
        leaves_[leafI]->slotPtr()->containedTriangles_;
    constRow elmts = containedTriangles[leaves_[leafI]->containedElements()];
    forAll(elmts, elI)
        triangles.append(elmts[elI]);
}

inline bool meshOctree::hasContainedEdges(const label leafI) const
{
    if( leaves_[leafI]->hasContainedEdges() )
        return true;

    return false;
}

inline void meshOctree::containedEdges
(
    const label leafI,
    DynList<label>& edges
) const
{
    edges.clear();

    if( !leaves_[leafI]->hasContainedEdges() )
    {
        edges.clear();
        return;
    }

    const VRWGraph& containedEdges =
        leaves_[leafI]->slotPtr()->containedEdges_;
    constRow elmts = containedEdges[leaves_[leafI]->containedEdges()];
    forAll(elmts, eI)
        edges.append(elmts[eI]);
}

inline const triSurf& meshOctree::surface() const
{
    return surface_;
}

inline const boundBox& meshOctree::rootBox() const
{
    return rootBox_;
}

inline const FixedList<meshOctreeCubeCoordinates, 26>&
meshOctree::regularityPositions() const
{
    return regularityPositions_;
}

inline const FixedList<FixedList<meshOctreeCubeCoordinates, 8>, 8>&
meshOctree::positionsOfLeavesAtNodes() const
{
    return vrtLeavesPos_;
}

inline const labelList& meshOctree::neiProcs() const
{
    return neiProcs_;
}

inline label meshOctree::findNeighbourOverNode
(
    const label leafI,
    const label nodeI
) const
{
    return findNeighbourOverNode(leaves_[leafI]->coordinates(), nodeI);
}

inline void meshOctree::findNeighboursOverEdge
(
    const label leafI,
    const label eI,
    DynList<label>& neighbourLeaves
) const
{
    findNeighboursOverEdge(leaves_[leafI]->coordinates(), eI, neighbourLeaves);
}

inline void meshOctree::findNeighboursInDirection
(
    const label leafI,
    const label dir,
    DynList<label>& neiLeaves
) const
{
    findNeighboursInDirection(leaves_[leafI]->coordinates(), dir, neiLeaves);
}

inline void meshOctree::findNeighboursForLeaf
(
    const label leafI,
    DynList<label>& neighbourLeaves
) const
{
    findNeighboursForLeaf(leaves_[leafI]->coordinates(), neighbourLeaves);
}

inline void meshOctree::findAllLeafNeighbours
(
    const label leafI,
    DynList<label>& neighbourLeaves
) const
{
    findAllLeafNeighbours(leaves_[leafI]->coordinates(), neighbourLeaves);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
