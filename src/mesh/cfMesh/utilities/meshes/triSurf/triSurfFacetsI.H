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

#include "triSurfFacets.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline label triSurfFacets::size() const
{
    return triangles_.size();
}

inline const LongList<labelledTri>& triSurfFacets::facets() const
{
    return triangles_;
}

inline const geometricSurfacePatchList& triSurfFacets::patches() const
{
    return patches_;
}

inline void triSurfFacets::appendTriangle(const labelledTri& tri)
{
    triangles_.append(tri);
}

inline void triSurfFacets::addFacetToSubset
(
    const label setI,
    const label triI
)
{
    Map<meshSubset>::iterator it = facetSubsets_.find(setI);
    if( it == facetSubsets_.end() )
        return;

    it().addElement(triI);
}

inline void triSurfFacets::removeFacetFromSubset
(
    const label setI,
    const label triI
)
{
    Map<meshSubset>::iterator it = facetSubsets_.find(setI);
    if( it == facetSubsets_.end() )
        return;

    it().removeElement(triI);
}

inline void triSurfFacets::facetInSubsets
(
    const label triI,
    DynList<label>& facetSubsets
) const
{
    facetSubsets.clear();

    forAllConstIter
    (
        Map<meshSubset>,
        facetSubsets_,
        it
    )
    {
        if( it().contains(triI) )
            facetSubsets.append(it.key());
    }
}

inline void triSurfFacets::facetSubsetIndices(DynList<label>& indices) const
{
    indices.clear();

    forAllConstIter
    (
        Map<meshSubset>,
        facetSubsets_,
        it
    )
        indices.append(it.key());
}

template<class ListType>
inline void triSurfFacets::facetsInSubset
(
    const label setI,
    ListType& facetsLabels
) const
{
    facetsLabels.clear();

    Map<meshSubset>::const_iterator it = facetSubsets_.find(setI);
    if( it == facetSubsets_.end() )
        return;

    it().containedElements(facetsLabels);
}

template<class ListType>
inline void triSurfFacets::updateFacetsSubsets(const ListType& newFacetsLabels)
{
    for
    (
        Map<meshSubset>::iterator it=facetSubsets_.begin();
        it!=facetSubsets_.end();
        ++it
    )
        it().updateSubset(newFacetsLabels);
}

inline const labelledTri& triSurfFacets::operator[](const label triI) const
{
    return triangles_[triI];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
