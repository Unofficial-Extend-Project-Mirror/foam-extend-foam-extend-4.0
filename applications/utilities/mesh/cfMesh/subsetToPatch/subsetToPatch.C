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
    Creates surface patches from surface subsets

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "triSurf.H"
#include "demandDrivenData.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void makePatchFromSubset
(
    triSurf& origSurf,
    const DynList<word>& subsetNames
)
{
    //- create new list of patches
    geometricSurfacePatchList newPatches
    (
        origSurf.patches().size() + subsetNames.size()
    );

    //- set names of the new patches
    forAll(origSurf.patches(), patchI)
        newPatches[patchI].name() = origSurf.patches()[patchI].name();

    forAll(subsetNames, subsetI)
        newPatches[origSurf.patches().size()+subsetI].name() =
            subsetNames[subsetI];

    //- create new triangles
    LongList<labelledTri> newTriangles(origSurf.facets());

    //- set patches for all triangles
    forAll(subsetNames, subsetI)
    {
        const label subsetID = origSurf.facetSubsetIndex(subsetNames[subsetI]);

        labelLongList subsetFaces;
        origSurf.facetsInSubset(subsetID, subsetFaces);

        const label regionI = origSurf.patches().size() + subsetI;

        forAll(subsetFaces, fI)
        {
            newTriangles[subsetFaces[fI]].region() = regionI;
        }
    }

    //- remove patches with no elements
    labelList nTrianglesInPatch(newPatches.size(), 0);
    forAll(newTriangles, triI)
        ++nTrianglesInPatch[newTriangles[triI].region()];

    Map<label> newPatchLabel;
    label counter(0);
    forAll(nTrianglesInPatch, patchI)
    {
        if( nTrianglesInPatch[patchI] )
            newPatchLabel.insert(patchI, counter++);
    }

    geometricSurfacePatchList copyPatches(counter);
    counter = 0;
    forAll(newPatches, patchI)
    {
        if( newPatchLabel.found(patchI) )
        {
            copyPatches[newPatchLabel[patchI]].name() =
                newPatches[patchI].name();
        }
    }

    newPatches = copyPatches;

    //- renumber the patches in the list of triangles
    forAll(newTriangles, triI)
        newTriangles[triI].region() =
            newPatchLabel[newTriangles[triI].region()];

    //- delete subsets converted to patches
    forAll(subsetNames, subsetI)
    {
        const label subsetID = origSurf.facetSubsetIndex(subsetNames[subsetI]);
        origSurf.removeFacetSubset(subsetID);
    }

    //- update subsets
    origSurf.updateFacetsSubsets(newPatchLabel);
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("input surface file");
    argList::validArgs.append("subset");
    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    word subsetName(args.args()[2]);

    triSurf* origSurfPtr = new triSurf(inFileName);

    DynList<word> subsetNames;
    const label subsetID = origSurfPtr->facetSubsetIndex(subsetName);
    if( subsetID >= 0 )
    {
        Warning << "Subset " << subsetName
        << " checking subsets containing this string!" << endl;

        DynList<label> existingSubsets;
        origSurfPtr->facetSubsetIndices(existingSubsets);

        forAll(existingSubsets, subsetI)
        {
            const word sName =
                origSurfPtr->facetSubsetName(existingSubsets[subsetI]);

            if( sName.substr(0, subsetName.size()) == subsetName )
            {
                subsetNames.append(sName);
            }
        }

        Info << "Converting " << subsetNames.size() << " subsets" << endl;
    }
    else
    {
        subsetNames.append(subsetName);
    }

    makePatchFromSubset(*origSurfPtr, subsetNames);
    origSurfPtr->writeSurface(inFileName);
    deleteDemandDrivenData(origSurfPtr);

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
