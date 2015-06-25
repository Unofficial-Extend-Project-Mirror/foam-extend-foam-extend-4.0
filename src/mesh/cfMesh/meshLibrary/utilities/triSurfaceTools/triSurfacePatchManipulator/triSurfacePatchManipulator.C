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

#include "triSurfacePatchManipulator.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"
#include "checkMeshDict.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfacePatchManipulator::triSurfacePatchManipulator(const triSurf& surface)
:
    surf_(surface),
    featureEdges_(surf_.edges().size(), direction(0)),
    facetInPatch_(),
    nPatches_(),
    newPatchNames_(),
    newPatchTypes_()
{
    allocateFeatureEdges();

    createPatches();
}

triSurfacePatchManipulator::~triSurfacePatchManipulator()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfacePatchManipulator::detectedSurfaceRegions
(
    VRWGraph& graph
) const
{
    graph.setSize(nPatches_);

    labelLongList nFacetsInPatch(nPatches_, 0);

    forAll(facetInPatch_, triI)
        ++nFacetsInPatch[facetInPatch_[triI]];

    graph.setSizeAndRowSize(nFacetsInPatch);

    nFacetsInPatch = 0;
    forAll(facetInPatch_, triI)
    {
        const label patchI = facetInPatch_[triI];

        graph(patchI, nFacetsInPatch[patchI]) = triI;
        ++nFacetsInPatch[patchI];
    }
}

const triSurf* triSurfacePatchManipulator::surfaceWithPatches
(
    IOdictionary* meshDictPtr,
    const word prefix,
    const bool forceOverwrite
) const
{
    //- collect patch information
    VRWGraph facetsInPatch;
    detectedSurfaceRegions(facetsInPatch);

    //- create new list of boundary patches
    LongList<labelledTri> newTriangles(facetInPatch_.size());
    label counter(0);
    geometricSurfacePatchList newPatches(nPatches_);

    if( forceOverwrite )
    {
        forAll(newPatches, patchI)
        {
            newPatches[patchI].name() = prefix+help::scalarToText(patchI);
            newPatches[patchI].geometricType() = "patch";
            newPatches[patchI].index() = patchI;
        }
    }
    else
    {
        forAll(facetsInPatch, patchI)
        {
            forAllRow(facetsInPatch, patchI, fpI)
            {
                const label origPatchI =
                    surf_[facetsInPatch(patchI, fpI)].region();
                newPatches[patchI].name() =
                    surf_.patches()[origPatchI].name() + '_' +
                    help::scalarToText(patchI);
                newPatches[patchI].geometricType() =
                    surf_.patches()[origPatchI].geometricType();
                newPatches[patchI].index() = patchI;
            }
        }
    }

    //- create triangles for the new surface
    labelLongList newFacetLabel(newTriangles.size(), -1);

    forAll(facetsInPatch, patchI)
        forAllRow(facetsInPatch, patchI, tI)
        {
            newFacetLabel[facetsInPatch(patchI, tI)] = counter;
            labelledTri tria = surf_[facetsInPatch(patchI, tI)];
            tria.region() = patchI;
            newTriangles[counter++] = tria;
        }

    //- create and return a new surface mesh
    triSurf* newSurfPtr =
        new triSurf
        (
            newTriangles,
            newPatches,
            edgeLongList(),
            surf_.points()
        );

    //- transfer facet subsets
    DynList<label> subsetIDs;
    surf_.facetSubsetIndices(subsetIDs);
    forAll(subsetIDs, subsetI)
    {
        const word sName = surf_.facetSubsetName(subsetIDs[subsetI]);

        const label newID = newSurfPtr->addFacetSubset(sName);

        labelLongList facetsInSubset;
        surf_.facetsInSubset(subsetIDs[subsetI], facetsInSubset);

        forAll(facetsInSubset, i)
        {
            const label fI = newFacetLabel[facetsInSubset[i]];

            newSurfPtr->addFacetToSubset(newID, fI);
        }
    }

    //- transfer point subsets
    surf_.pointSubsetIndices(subsetIDs);
    forAll(subsetIDs, subsetI)
    {
        const word sName = surf_.pointSubsetName(subsetIDs[subsetI]);

        const label newID = newSurfPtr->addPointSubset(sName);

        labelLongList pointsInSubset;
        surf_.pointsInSubset(subsetIDs[subsetI], pointsInSubset);

        forAll(pointsInSubset, i)
            newSurfPtr->addPointToSubset(newID, pointsInSubset[i]);
    }

    if( meshDictPtr )
    {
        //- create mapping between the patches on the original surface
        //- and the renamed patches
        std::map<word, wordList> patchesForPatch;
        std::map<word, word> patchTypes;

        const geometricSurfacePatchList& origPatches = surf_.patches();
        forAll(origPatches, patchI)
            patchTypes[origPatches[patchI].name()] =
                origPatches[patchI].geometricType();

        //- find the mapping of patch ids
        List<labelHashSet> patchToNewPatches(origPatches.size());
        forAll(facetsInPatch, patchI)
        {
            forAllRow(facetsInPatch, patchI, fI)
            {
                const label opatchI = surf_[facetsInPatch(patchI, fI)].region();

                patchToNewPatches[opatchI].insert(patchI);
            }
        }

        forAll(patchToNewPatches, patchI)
        {
            const word& pName = origPatches[patchI].name();
            patchesForPatch[pName].setSize
            (
                patchToNewPatches[patchI].size()
            );

            label counter(0);

            forAllConstIter(labelHashSet, patchToNewPatches[patchI], it)
                patchesForPatch[pName][counter++] = newPatches[it.key()].name();
        }

        //- update the values in meshDict based on the created patches
        checkMeshDict(*meshDictPtr).updateDictionaries
        (
            patchesForPatch,
            patchTypes
        );
    }

    return newSurfPtr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
