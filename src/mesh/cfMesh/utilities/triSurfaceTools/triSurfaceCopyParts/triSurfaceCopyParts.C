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

    //- mark patches which will be removed
    boolList removePatch(patches.size(), false);

    forAll(patches, patchI)
    {
        const word name = patches[patchI].name();

        forAll(parts, partI)
        {
            if( parts[partI] == name )
            {
                removePatch[patchI] = true;
                break;
            }
        }
    }

    //- select facets affected by the deletion of a patch
    forAll(surf_, triI)
    {
        if( removePatch[surf_[triI].region()] )
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
    LongList<labelledTri>& newTriangles = sm.facetsAccess();
    newTriangles.clear();
    sm.featureEdgesAccess().clear();
    sm.patchesAccess().setSize(1);
    sm.patchesAccess()[0] = geometricSurfacePatch("patch0", "patch", 0);

    //- copy selected patches
    forAll(copyFacets, triI)
    {
        if( !copyFacets[triI] )
            continue;

        const labelledTri& tri = surf_[triI];

        labelledTri newTri;
        newTri.region() = 0;

        forAll(tri, pI)
        {
            if( newPointLabel[tri[pI]] == -1 )
            {
                newPointLabel[tri[pI]] = nPoints;
                ++nPoints;
            }

            newTri[pI] = newPointLabel[tri[pI]];
        }

        newTriangles.append(newTri);
    }

    Info << "Copied triangles " << newTriangles.size() << endl;
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
