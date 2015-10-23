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

#include "triSurfaceCopyParts.H"
#include "triSurfModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceCopyParts::markFacetsForCopying
(
    const wordList& parts,
    boolList& copyFacets
) const
{
    copyFacets.setSize(surf_.size());
    copyFacets = false;

    const geometricSurfacePatchList& patches = surf_.patches();

    //- mark patches which will be copied
    boolList copyPatch(patches.size(), false);

    forAll(patches, patchI)
    {
        const word name = patches[patchI].name();

        forAll(parts, partI)
        {
            if( parts[partI] == name )
            {
                copyPatch[patchI] = true;
                break;
            }
        }
    }

    //- select facets affected by the deletion of a patch
    forAll(surf_, triI)
    {
        if( copyPatch[surf_[triI].region()] )
            copyFacets[triI] = true;
    }

    //- mark facets contained in selected subsets
    DynList<label> facetSubsetsIDs;
    surf_.facetSubsetIndices(facetSubsetsIDs);

    forAll(facetSubsetsIDs, i)
    {
        const word fsName = surf_.facetSubsetName(facetSubsetsIDs[i]);

        forAll(parts, partI)
        {
            if( parts[partI] == fsName )
            {
                labelLongList containedFacets;
                surf_.facetsInSubset(facetSubsetsIDs[i], containedFacets);

                forAll(containedFacets, cfI)
                    copyFacets[containedFacets[cfI]] = true;

                break;
            }
        }
    }
}

void triSurfaceCopyParts::copySurfaceMesh
(
    const boolList& copyFacets,
    triSurf& s
) const
{
    Info << "Starting copying surface parts" << endl;

    const pointField& pts = surf_.points();

    labelLongList newPointLabel(pts.size(), -1);

    label nPoints(0);

    //- create the modifier and delete data if there is any
    triSurfModifier sm(s);
    sm.facetsAccess().clear();
    sm.featureEdgesAccess().clear();
    sm.patchesAccess().setSize(surf_.patches().size());
    forAll(surf_.patches(), patchI)
        sm.patchesAccess()[patchI] = surf_.patches()[patchI];

    //- copy selected patches
    labelLongList newTriangleLabel(surf_.size(), -1);
    forAll(copyFacets, triI)
    {
        if( !copyFacets[triI] )
            continue;

        const labelledTri& tri = surf_[triI];

        labelledTri newTri;
        newTri.region() = tri.region();

        forAll(tri, pI)
        {
            if( newPointLabel[tri[pI]] == -1 )
            {
                newPointLabel[tri[pI]] = nPoints;
                ++nPoints;
            }

            newTri[pI] = newPointLabel[tri[pI]];
        }

        newTriangleLabel[triI] = s.size();
        s.appendTriangle(newTri);
    }

    Info << "Copied triangles " << s.size() << endl;
    Info << "Number of vertices " << nPoints << endl;

    //- copy vertices
    pointField& newPts = sm.pointsAccess();
    newPts.setSize(nPoints);

    forAll(newPointLabel, i)
    {
        if( newPointLabel[i] < 0 )
            continue;

        newPts[newPointLabel[i]] = pts[i];
    }

    //- copy point subsets
    DynList<label> subsetIds;
    surf_.pointSubsetIndices(subsetIds);
    forAll(subsetIds, subsetI)
    {
        const label origId = subsetIds[subsetI];
        const word sName = surf_.pointSubsetName(origId);

        labelLongList pointsInSubset;
        surf_.pointsInSubset(origId, pointsInSubset);

        const label newId = s.addPointSubset(sName);
        forAll(pointsInSubset, i)
        {
            const label newPointI = newPointLabel[pointsInSubset[i]];

            if( newPointI < 0 )
                continue;

            s.addPointToSubset(newId, newPointI);
        }
    }

    //- copy facet subsets
    surf_.facetSubsetIndices(subsetIds);
    forAll(subsetIds, subsetI)
    {
        const label origId = subsetIds[subsetI];
        const word sName = surf_.facetSubsetName(origId);

        labelLongList trianglesInSubset;
        surf_.facetsInSubset(origId, trianglesInSubset);

        const label newId = s.addFacetSubset(sName);
        forAll(trianglesInSubset, i)
        {
            const label newTriI = newTriangleLabel[trianglesInSubset[i]];

            if( newTriI < 0 )
                continue;

            s.addFacetToSubset(newId, newTriI);
        }
    }

    //- copy feature edges
    labelLongList newEdgeLabel(surf_.nFeatureEdges(), -1);
    const VRWGraph& pointEdges = surf_.pointEdges();
    const edgeLongList& edges = surf_.edges();
    const VRWGraph& edgeFacets = surf_.edgeFacets();
    forAll(newEdgeLabel, edgeI)
    {
        const edge& e = surf_.featureEdges()[edgeI];
        label eI(-1);
        forAllRow(pointEdges, e.start(), peI)
        {
            const label eJ = pointEdges(e.start(), peI);
            if( edges[eJ] == e )
            {
                eI = eJ;
                break;
            }
        }

        if( newPointLabel[e.start()] < 0 )
            continue;
        if( newPointLabel[e.end()] < 0 )
            continue;
        bool foundTriangle(false);
        forAllRow(edgeFacets, eI, efI)
        {
            if( newTriangleLabel[edgeFacets(eI, efI)] >= 0 )
            {
                foundTriangle = true;
                break;
            }
        }
        if( !foundTriangle )
            continue;

        newEdgeLabel[edgeI] = sm.featureEdgesAccess().size();
        sm.featureEdgesAccess().append
        (
            edge(newPointLabel[e.start()], newPointLabel[e.end()])
        );
    }

    //- copy subsets of feature edges
    surf_.edgeSubsetIndices(subsetIds);
    forAll(subsetIds, subsetI)
    {
        const label origId = subsetIds[subsetI];
        const word sName = surf_.edgeSubsetName(origId);

        labelLongList edgesInSubset;
        surf_.edgesInSubset(origId, edgesInSubset);

        const label newId = s.addEdgeSubset(sName);
        forAll(edgesInSubset, i)
        {
            const label newEdgeI = newEdgeLabel[edgesInSubset[i]];

            if( newEdgeI < 0 )
                continue;

            s.addEdgeToSubset(newId, newEdgeI);
        }
    }

    Info << "Finished copying surface parts" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfaceCopyParts::triSurfaceCopyParts(const triSurf& surface)
:
    surf_(surface)
{}

triSurfaceCopyParts::~triSurfaceCopyParts()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceCopyParts::copySurface(const wordList& patches, triSurf& s) const
{
    boolList copyFacets(surf_.size(), false);

    markFacetsForCopying(patches, copyFacets);

    copySurfaceMesh(copyFacets, s);
}

triSurf* triSurfaceCopyParts::copySurface(const wordList& patches) const
{
    boolList copyFacets(surf_.size(), false);

    markFacetsForCopying(patches, copyFacets);

    triSurf* sPtr = new triSurf();

    copySurfaceMesh(copyFacets, *sPtr);

    return sPtr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
