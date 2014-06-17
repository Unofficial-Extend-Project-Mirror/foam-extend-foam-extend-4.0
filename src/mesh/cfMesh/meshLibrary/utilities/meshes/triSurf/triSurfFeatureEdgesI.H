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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triSurfFeatureEdges.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline label triSurfFeatureEdges::nFeatureEdges() const
{
    return featureEdges_.size();
}

inline const edgeLongList& triSurfFeatureEdges::featureEdges() const
{
    return featureEdges_;
}

inline void triSurfFeatureEdges::appendFeatureEdge(const edge& e)
{
    featureEdges_.append(e);
}

inline void triSurfFeatureEdges::addEdgeToSubset
(
    const label setI,
    const label eI
)
{
    Map<meshSubset>::iterator it = featureEdgeSubsets_.find(setI);
    if( it == featureEdgeSubsets_.end() )
        return;

    it().addElement(eI);
}

inline void triSurfFeatureEdges::removeEdgeFromSubset
(
    const label setI,
    const label eI
)
{
    Map<meshSubset>::iterator it = featureEdgeSubsets_.find(setI);
    if( it == featureEdgeSubsets_.end() )
        return;

    it().removeElement(eI);
}

inline void triSurfFeatureEdges::edgeInSubsets
(
    const label eI,
    DynList<label>& edgeSubsets
) const
{
    edgeSubsets.clear();

    forAllConstIter
    (
        Map<meshSubset>,
        featureEdgeSubsets_,
        it
    )
    {
        if( it().contains(eI) )
            edgeSubsets.append(it.key());
    }
}

inline void triSurfFeatureEdges::edgeSubsetIndices
(
    DynList<label>& indices
) const
{
    indices.clear();

    forAllConstIter
    (
        Map<meshSubset>,
        featureEdgeSubsets_,
        it
    )
        indices.append(it.key());
}

template<class ListType>
inline void triSurfFeatureEdges::edgesInSubset
(
    const label setI,
    ListType& edgeLabels
) const
{
    edgeLabels.clear();

    Map<meshSubset>::const_iterator it = featureEdgeSubsets_.find(setI);
    if( it == featureEdgeSubsets_.end() )
        return;

    it().containedElements(edgeLabels);
}

template<class ListType>
inline void triSurfFeatureEdges::updateEdgeSubsets
(
    const ListType& newEdgeLabels
)
{
    for
    (
        Map<meshSubset>::iterator it=featureEdgeSubsets_.begin();
        it!=featureEdgeSubsets_.end();
        ++it
    )
        it().updateSubset(newEdgeLabels);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
