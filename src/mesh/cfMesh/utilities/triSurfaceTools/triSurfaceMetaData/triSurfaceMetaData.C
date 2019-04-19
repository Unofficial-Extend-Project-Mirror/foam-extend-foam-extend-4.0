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

#include "triSurfaceMetaData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceMetaData::createMetaData()
{
    metaDict_.clear();

    metaDict_.add("nPoints", surf_.points().size());
    metaDict_.add("nFacets", surf_.facets().size());
    metaDict_.add("nPatches", surf_.patches().size());
    metaDict_.add("nFeatureEdges", surf_.featureEdges().size());

    dictionary dict;

    //- store nformation about surface patches
    labelList nInPatch(surf_.patches().size(), 0);
    forAll(surf_, triI)
        ++nInPatch[surf_[triI].region()];

    forAll(surf_.patches(), patchI)
    {
        const geometricSurfacePatch& patch = surf_.patches()[patchI];

        dictionary pDict;
        pDict.add("type", patch.geometricType());
        pDict.add("nFacets", nInPatch[patchI]);

        dict.add(patch.name(), pDict);
    }

    metaDict_.add("patches", dict);

    //- store information about point subsets
    DynList<label> subsetIds;
    surf_.pointSubsetIndices(subsetIds);
    dict.clear();
    forAll(subsetIds, i)
    {
        dictionary sDict;

        labelLongList inSubset;
        surf_.pointsInSubset(subsetIds[i], inSubset);
        sDict.add("nPoints", inSubset.size());

        dict.add(surf_.pointSubsetName(subsetIds[i]), sDict);
    }

    metaDict_.add("pointSubsets", dict);

    //- store information about facet subsets
    subsetIds.clear();
    surf_.facetSubsetIndices(subsetIds);
    dict.clear();
    forAll(subsetIds, i)
    {
        dictionary sDict;

        labelLongList inSubset;
        surf_.facetsInSubset(subsetIds[i], inSubset);
        sDict.add("nFacets", inSubset.size());

        dict.add(surf_.facetSubsetName(subsetIds[i]), sDict);
    }

    metaDict_.add("facetSubsets", dict);

    //- store information about feature edge subsets
    subsetIds.clear();
    surf_.edgeSubsetIndices(subsetIds);
    dict.clear();
    forAll(subsetIds, i)
    {
        dictionary sDict;

        labelLongList inSubset;
        surf_.edgesInSubset(subsetIds[i], inSubset);
        sDict.add("nEdges", inSubset.size());

        dict.add(surf_.edgeSubsetName(subsetIds[i]), sDict);
    }

    metaDict_.add("featureEdgeSubsets", dict);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfaceMetaData::triSurfaceMetaData(const triSurf& surface)
:
    surf_(surface),
    metaDict_()
{
    createMetaData();
}

triSurfaceMetaData::~triSurfaceMetaData()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
