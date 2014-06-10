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

#include "triSurfaceRemoveFacets.H"
#include "triSurfModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceRemoveFacets::markFacetsForRemoval(boolList& removeFacet) const
{
    removeFacet.setSize(surf_.size());
    removeFacet = false;

    const geometricSurfacePatchList& patches = surf_.patches();

    //- mark patches which will be removed
    boolList removePatch(patches.size(), false);

    forAll(patches, patchI)
    {
        if( selectedEntities_.contains(patches[patchI].name()) )
            removePatch[patchI] = true;
    }

    //- select facets affected by the deletion of a patch
    forAll(surf_, triI)
    {
        if( removePatch[surf_[triI].region()] )
            removeFacet[triI] = true;
    }

    //- mark facets contained in selected subsets
    DynList<label> facetSubsetsIDs;
    surf_.facetSubsetIndices(facetSubsetsIDs);

    forAll(facetSubsetsIDs, i)
    {
        const word fsName = surf_.facetSubsetName(facetSubsetsIDs[i]);

        if( selectedEntities_.contains(fsName) )
        {
            labelLongList containedFacets;
            surf_.facetsInSubset(facetSubsetsIDs[i], containedFacets);

            forAll(containedFacets, cfI)
                removeFacet[containedFacets[cfI]] = true;
        }
    }
}

void triSurfaceRemoveFacets::removeFacets()
{
    boolList removeFacet;
    markFacetsForRemoval(removeFacet);

    //- calculate new indices of vertices and facets
    const pointField& points = surf_.points();
    labelLongList newPointLabel(surf_.points().size(), -1);
    labelLongList newFacetLabel(surf_.size(), -1);

    label pointCounter(0), facetCounter(0);

    forAll(removeFacet, triI)
    {
        if( removeFacet[triI] )
            continue;

        const labelledTri& tri = surf_[triI];

        forAll(tri, pI)
        {
            if( newPointLabel[tri[pI]] == -1 )
                newPointLabel[tri[pI]] = pointCounter++;
        }

        newFacetLabel[triI] = facetCounter++;
    }

    //- remove vertices
    pointField newPts(pointCounter);
    forAll(newPointLabel, pI)
    {
        if( newPointLabel[pI] < 0 )
            continue;

        newPts[newPointLabel[pI]] = points[pI];
    }

    triSurfModifier(surf_).pointsAccess().transfer(newPts);
    surf_.updatePointSubsets(newPointLabel);

    //- remove facets
    LongList<labelledTri> newFacets(facetCounter);

    forAll(newFacetLabel, triI)
    {
        if( newFacetLabel[triI] < 0 )
            continue;

        const labelledTri& tri = surf_[triI];

        newFacets[newFacetLabel[triI]] =
            labelledTri
            (
                newPointLabel[tri[0]],
                newPointLabel[tri[1]],
                newPointLabel[tri[2]],
                tri.region()
            );
    }

    triSurfModifier(surf_).facetsAccess().transfer(newFacets);
    surf_.updateFacetsSubsets(newFacetLabel);

    //- update feature edges
    const edgeLongList& featureEdges = surf_.featureEdges();
    label edgeCounter(0);
    labelLongList newFeatureEdgeLabel(featureEdges.size(), -1);

    forAll(featureEdges, feI)
    {
        const edge& e = featureEdges[feI];

        if( (newPointLabel[e.start()] < 0) || (newPointLabel[e.end()] < 0) )
            continue;

        newFeatureEdgeLabel[feI] = edgeCounter++;
    }

    edgeLongList newFeatureEdges(edgeCounter);
    forAll(newFeatureEdgeLabel, eI)
    {
        if( newFeatureEdgeLabel[eI] < 0 )
            continue;

        const edge& e = featureEdges[eI];

        newFeatureEdges[newFeatureEdgeLabel[eI]] =
            edge
            (
                newPointLabel[e.start()],
                newPointLabel[e.end()]
            );
    }

    triSurfModifier(surf_).featureEdgesAccess().transfer(newFeatureEdges);
    surf_.updateEdgeSubsets(newFeatureEdgeLabel);

    selectedEntities_.clear();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
