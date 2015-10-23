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

\*---------------------------------------------------------------------------*/

#include "surfaceMeshGeometryModification.H"
#include "demandDrivenData.H"
#include "dictionary.H"
#include "triSurf.H"

namespace Foam
{

// * * * * * * * * * * * * * * Private member functions* * * * * * * * * * * //

void surfaceMeshGeometryModification::checkModification()
{
    if( meshDict_.found("anisotropicSources") )
    {
        modificationActive_ = true;

        const dictionary& anisotropicDict =
            meshDict_.subDict("anisotropicSources");

        coordinateModifierPtr_ = new coordinateModifier(anisotropicDict);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceMeshGeometryModification::surfaceMeshGeometryModification
(
    const triSurf& surf,
    const dictionary& meshDict
)
:
    surf_(surf),
    meshDict_(meshDict),
    coordinateModifierPtr_(NULL),
    modificationActive_(false)
{
    checkModification();
}

surfaceMeshGeometryModification::~surfaceMeshGeometryModification()
{
    deleteDemandDrivenData(coordinateModifierPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool surfaceMeshGeometryModification::activeModification() const
{
    return modificationActive_;
}

const triSurf* surfaceMeshGeometryModification::modifyGeometry() const
{
    if( !modificationActive_ )
    {
        WarningIn
        (
            "const triSurf* surfaceMeshGeometryModification"
            "::modifyGeometry() const"
        ) << "Modification is not active" << endl;

        return NULL;
    }

    const pointField& pts = surf_.points();

    pointField newPts(pts.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(pts, pointI)
        newPts[pointI] = coordinateModifierPtr_->modifiedPoint(pts[pointI]);

    triSurf* newSurf =
        new triSurf
        (
            surf_.facets(),
            surf_.patches(),
            surf_.featureEdges(),
            newPts
        );

    //- copy subsets
    DynList<label> sIds;

    //- copy facet subsets
    surf_.facetSubsetIndices(sIds);
    forAll(sIds, i)
    {
        const label newId =
            newSurf->addFacetSubset(surf_.facetSubsetName(sIds[i]));

        labelLongList facetsInSubset;
        surf_.facetsInSubset(sIds[i], facetsInSubset);

        forAll(facetsInSubset, fI)
            newSurf->addFacetToSubset(newId, facetsInSubset[fI]);
    }

    //- copy point subsets
    surf_.pointSubsetIndices(sIds);
    forAll(sIds, i)
    {
        const label newId =
            newSurf->addPointSubset(surf_.pointSubsetName(sIds[i]));

        labelLongList pointsInSubset;
        surf_.pointsInSubset(sIds[i], pointsInSubset);

        forAll(pointsInSubset, pI)
            newSurf->addPointToSubset(newId, pointsInSubset[pI]);
    }

    //- copy subsets of feature edges
    surf_.edgeSubsetIndices(sIds);
    forAll(sIds, i)
    {
        const label newId =
            newSurf->addEdgeSubset(surf_.edgeSubsetName(sIds[i]));

        labelLongList edgesInSubset;
        surf_.edgesInSubset(sIds[i], edgesInSubset);

        forAll(edgesInSubset, eI)
            newSurf->addEdgeToSubset(newId, edgesInSubset[eI]);
    }

    return newSurf;
}

const triSurf* surfaceMeshGeometryModification::
revertGeometryModification() const
{
    if( !modificationActive_ )
    {
        WarningIn
        (
            "const triSurf* surfaceMeshGeometryModification"
            "::revertGeometryModification() const"
        ) << "Modification is not active" << endl;

        return NULL;
    }

    const pointField& pts = surf_.points();

    pointField newPts(pts.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(pts, pointI)
        newPts[pointI] =
            coordinateModifierPtr_->backwardModifiedPoint(pts[pointI]);

    triSurf* newSurf =
        new triSurf
        (
            surf_.facets(),
            surf_.patches(),
            surf_.featureEdges(),
            newPts
        );

    //- copy subsets
    DynList<label> sIds;

    //- copy facet subsets
    surf_.facetSubsetIndices(sIds);
    forAll(sIds, i)
    {
        const label newId =
            newSurf->addFacetSubset(surf_.facetSubsetName(sIds[i]));

        labelLongList facetsInSubset;
        surf_.facetsInSubset(sIds[i], facetsInSubset);

        forAll(facetsInSubset, fI)
            newSurf->addFacetToSubset(newId, facetsInSubset[fI]);
    }

    //- copy point subsets
    surf_.pointSubsetIndices(sIds);
    forAll(sIds, i)
    {
        const label newId =
            newSurf->addPointSubset(surf_.pointSubsetName(sIds[i]));

        labelLongList pointsInSubset;
        surf_.pointsInSubset(sIds[i], pointsInSubset);

        forAll(pointsInSubset, pI)
            newSurf->addPointToSubset(newId, pointsInSubset[pI]);
    }

    //- copy subsets of feature edges
    surf_.edgeSubsetIndices(sIds);
    forAll(sIds, i)
    {
        const label newId =
            newSurf->addEdgeSubset(surf_.edgeSubsetName(sIds[i]));

        labelLongList edgesInSubset;
        surf_.edgesInSubset(sIds[i], edgesInSubset);

        forAll(edgesInSubset, eI)
            newSurf->addEdgeToSubset(newId, edgesInSubset[eI]);
    }

    return newSurf;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
