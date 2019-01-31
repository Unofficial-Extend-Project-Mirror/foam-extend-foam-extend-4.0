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

#include "checkMeshDict.H"
#include "patchRefinementList.H"
#include "PtrList.H"
#include "LongList.H"
#include "objectRefinement.H"
#include "coordinateModification.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkMeshDict::checkBasicSettings() const
{
    //- check if maxCellSize is valid
    const scalar maxCellSize = readScalar(meshDict_.lookup("maxCellSize"));

    if( maxCellSize < 0 )
        FatalErrorIn
        (
            "void checkMeshDict::checkBasicSettings() const"
        ) << "maxCellSize is negative! Cannot generate the mesh!!"
          << exit(FatalError);

    //- check if boundaryCellSize makes sense
    if( meshDict_.found("boundaryCellSize") )
    {
        const scalar bcs = readScalar(meshDict_.lookup("boundaryCellSize"));

        if( bcs < 0 )
        {
            WarningIn
            (
                "void checkMeshDict::checkBasicSettings() const"
            ) << "Boundary cell size is negative!!" << endl;
        }

        if( meshDict_.found("boundaryCellSizeRefinementThickness") )
        {
            const scalar thickness =
                readScalar
                (
                    meshDict_.lookup("boundaryCellSizeRefinementThickness")
                );

            if( thickness < 0 )
            {
                WarningIn
                (
                    "void checkMeshDict::checkBasicSettings() const"
                ) << "Boundary cell size refinement thickness is negative!!"
                  << endl;
            }
        }
    }

    //- check if minCellSize is valid
    if( meshDict_.found("minCellSize") )
    {
        const scalar mcs = readScalar(meshDict_.lookup("minCellSize"));

        if( mcs < 0 )
        {
            FatalErrorIn
            (
                "void checkMeshDict::checkBasicSettings() const"
            ) << "Minimum cell size for automatic refinement is negative!!"
              << exit(FatalError);
        }

    }

    //- check if keepCellsIntersectingBoundary can be read correctly
    if( meshDict_.found("keepCellsIntersectingBoundary") )
    {
        const bool keep =
            readBool(meshDict_.lookup("keepCellsIntersectingBoundary"));

        if( keep && meshDict_.found("checkForGluedMesh") )
        {
            readBool(meshDict_.lookup("checkForGluedMesh"));
        }
    }

    //- check if enforceConstraints is available
    if( meshDict_.found("enforceGeometryConstraints") )
    {
        readBool(meshDict_.lookup("enforceGeometryConstraints"));
    }
}

void checkMeshDict::checkPatchCellSize() const
{
    if( meshDict_.found("patchCellSize") )
    {
        if( meshDict_.isDict("patchCellSize") )
        {
            const dictionary& dict = meshDict_.subDict("patchCellSize");

            const wordList patchNames = dict.toc();
            patchNames.size();
        }
        else
        {
            patchRefinementList prl(meshDict_.lookup("patchCellSize"));
            prl.size();
        }
    }
}

void checkMeshDict::checkSubsetCellSize() const
{
    if( meshDict_.found("subsetCellSize") )
    {
        if( meshDict_.isDict("subsetCellSize") )
        {
            const dictionary& dict = meshDict_.subDict("subsetCellSize");

            const wordList subsetNames = dict.toc();
            subsetNames.size();
        }
        else
        {
            patchRefinementList prl(meshDict_.lookup("patchCellSize"));
        }
    }
}

void checkMeshDict::checkLocalRefinementLevel() const
{
    if( meshDict_.found("localRefinement") )
    {
        if( meshDict_.isDict("localRefinement") )
        {
            const dictionary& refDict = meshDict_.subDict("localRefinement");
            const wordList entries = refDict.toc();

            forAll(entries, dictI)
            {
                const dictionary& dict = refDict.subDict(entries[dictI]);

                if( dict.found("cellSize") )
                {
                    const scalar cs = readScalar(dict.lookup("cellSize"));

                    if( cs < 0.0 )
                    {
                        WarningIn
                        (
                        "void checkMeshDict::checkLocalRefinementLevel() const"
                        ) << "Cell size for " << entries[dictI]
                             << " is negative" << endl;
                    }
                }
                else if( dict.found("additionalRefinementLevels") )
                {
                    const label nLevels =
                        readLabel(dict.lookup("additionalRefinementLevels"));

                    if( nLevels < 0 )
                    {
                        WarningIn
                        (
                        "void checkMeshDict::checkLocalRefinementLevel() const"
                        ) << "Refinement level for " << entries[dictI]
                             << " is negative" << endl;
                    }
                }
                else
                {
                    FatalErrorIn
                    (
                        "void checkMeshDict::checkLocalRefinementLevel() const"
                    ) << "Cannot read keyword"
                      << " additionalRefinementLevels or cellSize"
                      << "for " << entries[dictI] << exit(FatalError);
                }

                if( dict.found("refinementThickness") )
                {
                    const scalar s =
                        readScalar(dict.lookup("refinementThickness"));

                    if( s < 0 )
                    {
                        WarningIn
                        (
                        "void checkMeshDict::checkLocalRefinementLevel() const"
                        ) << "Refinement thickness for " << entries[dictI]
                             << " is negative" << endl;
                    }
                }
            }
        }
        else
        {
            FatalErrorIn
            (
                "void checkMeshDict::checkLocalRefinementLevel() const"
            ) << "Cannot read localRefinement" << exit(FatalError);
        }
    }
}

void checkMeshDict::checkKeepCellsIntersectingPatches() const
{
    if( meshDict_.found("keepCellsIntersectingPatches") )
    {
        if( meshDict_.isDict("keepCellsIntersectingPatches") )
        {
            const dictionary& dict =
                meshDict_.subDict("keepCellsIntersectingPatches");

            const wordList patchNames = dict.toc();
            patchNames.size();
        }
        else
        {
            wordList kcip(meshDict_.lookup("keepCellsIntersectingPatches"));
        }
    }
}

void checkMeshDict::checkRemoveCellsIntersectingPatches() const
{
    if( meshDict_.found("removeCellsIntersectingPatches") )
    {
        if( meshDict_.isDict("removeCellsIntersectingPatches") )
        {
            const dictionary& dict =
                meshDict_.subDict("removeCellsIntersectingPatches");

            const wordList patchNames = dict.toc();
            patchNames.size();
        }
        else
        {
            wordList kcip(meshDict_.lookup("removeCellsIntersectingPatches"));
        }
    }
}

void checkMeshDict::checkObjectRefinements() const
{
    if( meshDict_.found("objectRefinements") )
    {
        PtrList<objectRefinement> refObjects;

        if( meshDict_.isDict("objectRefinements") )
        {
            const dictionary& dict = meshDict_.subDict("objectRefinements");
            const wordList objectNames = dict.toc();

            refObjects.setSize(objectNames.size());

            forAll(refObjects, objectI)
            {
                const entry& objectEntry =
                    dict.lookupEntry(objectNames[objectI], false, false);

                refObjects.set
                (
                    objectI,
                    objectRefinement::New
                    (
                        objectEntry.keyword(),
                        objectEntry.dict()
                    )
                );
            }
        }
        else
        {
            Istream& is = meshDict_.lookup("objectRefinements");

            PtrList<entry> objectEntries(is);
            refObjects.setSize(objectEntries.size());

            forAll(refObjects, objectI)
            {
                refObjects.set
                (
                    objectI,
                    objectRefinement::New
                    (
                        objectEntries[objectI].keyword(),
                        objectEntries[objectI].dict()
                    )
                );
            }
        }

        forAll(refObjects, oI)
        {
            if( refObjects[oI].refinementThickness() < 0.0 )
            {
                WarningIn
                (
                    "void checkMeshDict::checkObjectRefinements() const"
                ) << "Refinement thickness specified for object "
                  << refObjects[oI].name() << " is negative!!" << endl;
            }
        }
    }
}

void checkMeshDict::checkAnisotropicSources() const
{
    if( meshDict_.found("anisotropicSources") )
    {
        PtrList<coordinateModification> anisotropicObjects;

        if( meshDict_.isDict("anisotropicSources") )
        {
            const dictionary& dict = meshDict_.subDict("anisotropicSources");
            const wordList objectNames = dict.toc();

            anisotropicObjects.setSize(objectNames.size());

            forAll(anisotropicObjects, objectI)
            {
                const entry& objectEntry =
                    dict.lookupEntry(objectNames[objectI], false, false);

                anisotropicObjects.set
                (
                    objectI,
                    coordinateModification::New
                    (
                        objectEntry.keyword(),
                        objectEntry.dict()
                    )
                );
            }
        }
    }
}

void checkMeshDict::checkSurfaceRefinements() const
{
    if( meshDict_.found("surfaceMeshRefinement") )
    {
        const dictionary& surfaces = meshDict_.subDict("surfaceMeshRefinement");

        const wordList surfaceSources = surfaces.toc();

        forAll(surfaceSources, surfI)
        {
            if( surfaces.isDict(surfaceSources[surfI]) )
            {
                const dictionary& dict =
                    surfaces.subDict(surfaceSources[surfI]);

                if( dict.found("surfaceFile") )
                {
                    const fileName fName(dict.lookup("surfaceFile"));

                    if( !isFile(fName) )
                        FatalErrorIn
                        (
                        "void checkMeshDict::checkSurfaceRefinements() const"
                        ) << "Surface file " << fName
                          << " does not exist or is not readable!!"
                          << exit(FatalError);
                }
                else
                {
                    FatalErrorIn
                    (
                        "void checkMeshDict::checkSurfaceRefinements() const"
                    ) << "Missing surfaceFile for entry "
                      << surfaceSources[surfI] << exit(FatalError);
                }

                if( dict.found("cellSize") )
                {
                    const scalar cs = readScalar(dict.lookup("cellSize"));

                    if( cs < VSMALL )
                        FatalErrorIn
                        (
                            "void checkMeshDict::"
                            "checkSurfaceRefinements() const"
                        ) << "Cell size for entry " << surfaceSources[surfI]
                          << " is extremely small or negative!!"
                          << exit(FatalError);
                }
                else if( dict.found("additionalRefinementLevels") )
                {
                    const label nLev =
                        readLabel(dict.lookup("additionalRefinementLevels"));

                    if( nLev < 0 )
                    {
                        FatalErrorIn
                        (
                            "void checkMeshDict::"
                            "checkSurfaceRefinements() const"
                        ) << "Number refinement levels for entry "
                          << surfaceSources[surfI] << " is negative!!"
                          << exit(FatalError);
                    }
                }
                else
                {
                    FatalErrorIn
                    (
                        "void checkMeshDict::checkSurfaceRefinements() const"
                    ) << "Missing cellSize or additionalRefinementLevels"
                      << " for entry " << surfaceSources[surfI]
                      << exit(FatalError);
                }

                if( dict.found("refinementThickness") )
                {
                    const scalar cs =
                        readScalar(dict.lookup("refinementThickness"));

                    if( cs < VSMALL )
                        WarningIn
                        (
                            "void checkMeshDict::"
                            "checkSurfaceRefinements() const"
                        ) << "Refinement thickness for entry "
                          << surfaceSources[surfI]
                          << " is extremely small or negative!!" << endl;
                }
            }
            else
            {
                FatalErrorIn
                (
                    "void checkMeshDict::checkSurfaceRefinements() const"
                ) << "Dictionary " << surfaceSources[surfI]
                  << " does not exist!!"
                  << exit(FatalError);
            }
        }
    }
}

void checkMeshDict::checkEdgeMeshRefinements() const
{
    if( meshDict_.found("edgeMeshRefinement") )
    {
        const dictionary& edgeMeshes = meshDict_.subDict("edgeMeshRefinement");

        const wordList edgeMeshSources = edgeMeshes.toc();

        forAll(edgeMeshSources, emI)
        {
            if( edgeMeshes.isDict(edgeMeshSources[emI]) )
            {
                const dictionary& dict =
                    edgeMeshes.subDict(edgeMeshSources[emI]);

                if( dict.found("edgeFile") )
                {
                    const fileName fName(dict.lookup("edgeFile"));

                    if( !isFile(fName) )
                        FatalErrorIn
                        (
                        "void checkMeshDict::checkEdgeMeshRefinements() const"
                        ) << "Edge mesh file " << fName
                          << " does not exist or is not readable!!"
                          << exit(FatalError);
                }
                else
                {
                    FatalErrorIn
                    (
                        "void checkMeshDict::checkEdgeMeshRefinements() const"
                    ) << "Missing edgeFilw for entry "
                      << edgeMeshSources[emI] << exit(FatalError);
                }

                if( dict.found("cellSize") )
                {
                    const scalar cs = readScalar(dict.lookup("cellSize"));

                    if( cs < VSMALL )
                        FatalErrorIn
                        (
                            "void checkMeshDict::"
                            "checkEdgeMeshRefinements() const"
                        ) << "Cell size for entry " << edgeMeshSources[emI]
                          << " is extremely small or negative!!"
                          << exit(FatalError);
                }
                else if( dict.found("additionalRefinementLevels") )
                {
                    const label nLev =
                        readLabel(dict.lookup("additionalRefinementLevels"));

                    if( nLev < 0 )
                    {
                        FatalErrorIn
                        (
                            "void checkMeshDict::"
                            "checkEdgeMeshRefinements() const"
                        ) << "Number refinement levels for entry "
                          << edgeMeshSources[emI] << " is negative!!"
                          << exit(FatalError);
                    }
                }
                else
                {
                    FatalErrorIn
                    (
                        "void checkMeshDict::checkEdgeMeshRefinements() const"
                    ) << "Missing cellSize or additionalRefinementLevels"
                      << " for entry " << edgeMeshSources[emI]
                      << exit(FatalError);
                }

                if( dict.found("refinementThickness") )
                {
                    const scalar cs =
                        readScalar(dict.lookup("refinementThickness"));

                    if( cs < VSMALL )
                        WarningIn
                        (
                            "void checkMeshDict::"
                            "checkEdgeMeshRefinements() const"
                        ) << "Refinement thickness for entry "
                          << edgeMeshSources[emI]
                          << " is extremely small or negative!!" << endl;
                }
            }
            else
            {
                FatalErrorIn
                (
                    "void checkMeshDict::checkEdgeMeshRefinements() const"
                ) << "Dictionary " << edgeMeshSources[emI]
                  << " does not exist!!"
                  << exit(FatalError);
            }
        }
    }
}

void checkMeshDict::checkBoundaryLayers() const
{
    if( meshDict_.found("boundaryLayers") )
    {
        const dictionary& bndLayers = meshDict_.subDict("boundaryLayers");

        //- read global properties
        if( bndLayers.found("nLayers") )
        {
            readLabel(bndLayers.lookup("nLayers"));
        }
        if( bndLayers.found("thicknessRatio") )
        {
            readScalar(bndLayers.lookup("thicknessRatio"));
        }
        if( bndLayers.found("maxFirstLayerThickness") )
        {
            readScalar(bndLayers.lookup("maxFirstLayerThickness"));
        }

        //- patch-based properties
        if( bndLayers.isDict("patchBoundaryLayers") )
        {
            const dictionary& patchBndLayers =
                bndLayers.subDict("patchBoundaryLayers");
            const wordList patchNames = patchBndLayers.toc();

            forAll(patchNames, patchI)
            {
                const word pName = patchNames[patchI];

                if( patchBndLayers.isDict(pName) )
                {
                    const dictionary& patchDict =
                        patchBndLayers.subDict(pName);

                    if( patchDict.found("nLayers") )
                    {
                        readLabel(patchDict.lookup("nLayers"));
                    }
                    if( patchDict.found("thicknessRatio") )
                    {
                        readScalar(patchDict.lookup("thicknessRatio"));
                    }
                    if( patchDict.found("maxFirstLayerThickness") )
                    {
                        readScalar(patchDict.lookup("maxFirstLayerThickness"));
                    }
                    if( patchDict.found("allowDiscontinuity") )
                    {
                        readBool(patchDict.lookup("allowDiscontinuity"));
                    }
                }
                else
                {
                    Warning << "Cannot refine layer for patch "
                        << patchNames[patchI] << endl;
                }
            }
        }

        //- check for existence of boundary layer smoothing
        if( bndLayers.found("optimiseLayer") )
        {
            readBool(bndLayers.lookup("optimiseLayer"));
        }

        if( bndLayers.found("optimisationParameters") )
        {
            const dictionary& optParams =
                bndLayers.subDict("optimisationParameters");

            if( optParams.found("reCalculateNormals") )
            {
                readBool(optParams.lookup("reCalculateNormals"));
            }

            if( optParams.found("nSmoothNormals") )
            {
                const label nSmoothNormals =
                    readLabel(optParams.lookup("nSmoothNormals"));

                if( nSmoothNormals < 0 )
                    FatalErrorIn
                    (
                        "void checkMeshDict::checkBoundaryLayers() const"
                    ) << "nSmoothNormals must not be negative!"
                      << exit(FatalError);
            }

            if( optParams.found("featureSizeFactor") )
            {
                const scalar featureSizeFactor =
                    readScalar(optParams.lookup("featureSizeFactor"));

                if( featureSizeFactor >= 1.0 || featureSizeFactor < 0.0 )
                    FatalErrorIn
                    (
                        "void checkMeshDict::checkBoundaryLayers() const"
                    ) << "Feature size factor is out"
                      << " of a valid range 0 to 1" << exit(FatalError);
            }

            if( optParams.found("relThicknessTol") )
            {
                const scalar relThicknessTol =
                    readScalar(optParams.lookup("relThicknessTol"));

                if( relThicknessTol >= 1.0 || relThicknessTol < 0.0 )
                    FatalErrorIn
                    (
                        "void checkMeshDict::checkBoundaryLayers() const"
                    ) << "Relative thickness tolerance is out"
                      << " of a valid range 0 to 1" << exit(FatalError);
            }

            if( optParams.found("maxNumIterations") )
            {
                const label maxNumIterations =
                    readLabel(optParams.lookup("maxNumIterations"));

                if( maxNumIterations < 0 )
                    FatalErrorIn
                    (
                        "void checkMeshDict::checkBoundaryLayers() const"
                    ) << "maxNumIterations must not be negative!"
                      << exit(FatalError);
            }
        }
    }
}

void checkMeshDict::checkRenameBoundary() const
{
    if( meshDict_.found("renameBoundary") )
    {
        const dictionary& dict = meshDict_.subDict("renameBoundary");
        if( dict.found("newPatchNames") )
        {
            if( dict.isDict("newPatchNames") )
            {
                const dictionary& patchDicts = dict.subDict("newPatchNames");

                const wordList patchNames = patchDicts.toc();

                forAll(patchNames, patchI)
                {
                    const word& pName = patchNames[patchI];

                    if( !patchDicts.isDict(pName) )
                        FatalErrorIn
                        (
                            "void checkMeshDict::checkRenameBoundary() const"
                        ) << "Entry " << pName
                          << " is not a dictionary" << exit(FatalError);

                    const dictionary dict = patchDicts.subDict(pName);

                    if( !dict.found("newName") )
                        FatalErrorIn
                        (
                            "void checkMeshDict::checkRenameBoundary() const"
                        ) << "Dictionary " << pName
                          << " does not contain a newName keyword"
                          << exit(FatalError);
                }
            }
            else
            {
                const PtrList<entry> patchesToRename
                (
                    dict.lookup("newPatchNames")
                );

                forAll(patchesToRename, patchI)
                {
                    const word& pName = patchesToRename[patchI].keyword();

                    const dictionary dict = patchesToRename[patchI].dict();

                    if( !dict.found("newName") )
                        FatalErrorIn
                        (
                            "void checkMeshDict::checkRenameBoundary() const"
                        ) << "Dictionary " << pName
                          << " does not contain a newName keyword"
                          << exit(FatalError);
                }
            }
        }
    }
}

//Mesh quality criteria specified by the user
void checkMeshDict::checkQualitySettings() const
{

    if( meshDict_.found("meshQualitySettings") )
    {
        const dictionary& qualityDict = meshDict_.subDict("meshQualitySettings");

        //- read maximum non-orthogonality defined by the user
        if( qualityDict.found("maxNonOrthogonality") )
        {
            readScalar(qualityDict.lookup("maxNonOrthogonality"));
        }

        //- read maximum skewness defined by the user
        if( qualityDict.found("maxSkewness") )
        {
            readScalar(qualityDict.lookup("maxSkewness"));
        }

        //- read minimum volume of the face pyramid defined by the user
        if( qualityDict.found("minPyramidVolume") )
        {
            readScalar(qualityDict.lookup("minPyramidVolume"));
        }

        //- read face flatness defined by the user
        if( qualityDict.found("faceFlatness") )
        {
            readScalar(qualityDict.lookup("faceFlatness"));
        }

        //- read minimum tetrahedral part of a cell defined by the user
        if( qualityDict.found("minCellPartTetrahedra") )
        {
            readScalar(qualityDict.lookup("minCellPartTetrahedra"));
        }

        //- read minimum area of a face defined by the user
        if( qualityDict.found("minimumFaceArea") )
        {
            readScalar(qualityDict.lookup("minimumFaceArea"));
        }


    }
}

void checkMeshDict::checkEntries() const
{
    checkBasicSettings();

    checkPatchCellSize();

    checkSubsetCellSize();

    checkSurfaceRefinements();

    checkKeepCellsIntersectingPatches();

    checkRemoveCellsIntersectingPatches();

    checkObjectRefinements();

    checkAnisotropicSources();

    checkBoundaryLayers();

    checkRenameBoundary();

    checkQualitySettings();
}

void checkMeshDict::updatePatchCellSize
(
    const std::map<word, wordList>& patchesFromPatch
)
{
    if( meshDict_.found("patchCellSize") )
    {
        LongList<patchRefinement> updatedPatchRefinement;

        if( meshDict_.isDict("patchCellSize") )
        {
            const dictionary dict = meshDict_.subDict("patchCellSize");

            const wordList patchNames = dict.toc();

            forAll(patchNames, patchI)
            {
                const word& pName = patchNames[patchI];

                std::map<word, wordList>::const_iterator it =
                    patchesFromPatch.find(pName);
                if( it == patchesFromPatch.end() )
                    continue;

                const wordList& updatedPatchNames = it->second;

                const dictionary& pDict = dict.subDict(pName);
                const scalar cellSize = readScalar(pDict.lookup("cellSize"));

                forAll(updatedPatchNames, nameI)
                    updatedPatchRefinement.append
                    (
                        patchRefinement
                        (
                            updatedPatchNames[nameI],
                            cellSize
                        )
                    );
            }
        }
        else
        {
            patchRefinementList prl(meshDict_.lookup("patchCellSize"));
            forAll(prl, prlI)
            {
                const word& pName = prl[prlI].patchName();
                const scalar cellSize = prl[prlI].cellSize();

                std::map<word, wordList>::const_iterator it =
                    patchesFromPatch.find(pName);

                if( it == patchesFromPatch.end() )
                    continue;

                const wordList& updatedPatchNames = it->second;
                forAll(updatedPatchNames, nameI)
                    updatedPatchRefinement.append
                    (
                        patchRefinement
                        (
                            updatedPatchNames[nameI],
                            cellSize
                        )
                    );
            }
        }

        meshDict_.add("patchCellSize", updatedPatchRefinement, true);
    }
}

void checkMeshDict::updateSubsetCellSize
(
    const std::map<word, wordList>& /*patchesFromPatch*/
)
{

}

void checkMeshDict::updateLocalRefinement
(
    const std::map<word, wordList>& patchesFromPatch
)
{
    if( meshDict_.found("localRefinement") )
    {
        if( meshDict_.isDict("localRefinement") )
        {
            dictionary& dict = meshDict_.subDict("localRefinement");

            const wordList entries = dict.toc();

            forAll(entries, dictI)
            {
                const word& pName = entries[dictI];

                std::map<word, wordList>::const_iterator it =
                    patchesFromPatch.find(pName);
                if( it == patchesFromPatch.end() )
                    continue;

                const wordList& updatedPatchNames = it->second;

                const dictionary& pDict = dict.subDict(pName);
                dictionary copy = pDict;

                //- add new patches
                forAll(updatedPatchNames, nameI)
                    dict.add(updatedPatchNames[nameI], copy);

                //- remove the current patch
                dict.remove(pName);
            }
        }
    }
}

void checkMeshDict::updateKeepCellsIntersectingPatches
(
    const std::map<word, wordList>& patchesFromPatch
)
{
    if( meshDict_.found("keepCellsIntersectingPatches") )
    {
        LongList<word> updatedPatchNames;
        if( meshDict_.isDict("keepCellsIntersectingPatches") )
        {
            const dictionary& dict =
                meshDict_.subDict("keepCellsIntersectingPatches");

            const wordList patchNames = dict.toc();
            forAll(patchNames, patchI)
            {
                const word& pName = patchNames[patchI];

                std::map<word, wordList>::const_iterator it =
                    patchesFromPatch.find(pName);

                if( it == patchesFromPatch.end() )
                {
                    updatedPatchNames.append(pName);
                    continue;
                }

                const wordList& newPatchNames = it->second;

                forAll(newPatchNames, nameI)
                    updatedPatchNames.append(newPatchNames[nameI]);
            }
        }
        else
        {
            wordList kcip(meshDict_.lookup("keepCellsIntersectingPatches"));

            forAll(kcip, i)
            {
                const word& pName = kcip[i];

                std::map<word, wordList>::const_iterator it =
                    patchesFromPatch.find(pName);

                if( it == patchesFromPatch.end() )
                    updatedPatchNames.append(pName);

                const wordList& newPatchNames = it->second;

                forAll(newPatchNames, nameI)
                    updatedPatchNames.append(newPatchNames[nameI]);
            }
        }

        meshDict_.add("keepCellsIntersectingPatches", updatedPatchNames, true);
    }
}


void checkMeshDict::updateRemoveCellsIntersectingPatches
(
    const std::map<word, wordList>& patchesFromPatch
)
{
    if( meshDict_.found("removeCellsIntersectingPatches") )
    {
        LongList<word> updatedPatchNames;
        if( meshDict_.isDict("removeCellsIntersectingPatches") )
        {
            const dictionary& dict =
                meshDict_.subDict("removeCellsIntersectingPatches");

            const wordList patchNames = dict.toc();
            forAll(patchNames, patchI)
            {
                const word& pName = patchNames[patchI];

                std::map<word, wordList>::const_iterator it =
                    patchesFromPatch.find(pName);

                if( it == patchesFromPatch.end() )
                {
                    updatedPatchNames.append(pName);
                    continue;
                }

                const wordList& newPatchNames = it->second;

                forAll(newPatchNames, nameI)
                    updatedPatchNames.append(newPatchNames[nameI]);
            }
        }
        else
        {
            wordList kcip(meshDict_.lookup("removeCellsIntersectingPatches"));

            forAll(kcip, i)
            {
                const word& pName = kcip[i];

                std::map<word, wordList>::const_iterator it =
                    patchesFromPatch.find(pName);

                if( it == patchesFromPatch.end() )
                    updatedPatchNames.append(pName);

                const wordList& newPatchNames = it->second;

                forAll(newPatchNames, nameI)
                    updatedPatchNames.append(newPatchNames[nameI]);
            }
        }

        meshDict_.add
        (
            "removeCellsIntersectingPatches",
            updatedPatchNames,
            true
        );
    }
}

void checkMeshDict::updateObjectRefinements
(
    const std::map<word, wordList>& /*patchesFromPatch*/
)
{

}

void checkMeshDict::updateBoundaryLayers
(
    const std::map<word, wordList>& patchesFromPatch
)
{
    if( meshDict_.isDict("boundaryLayers") )
    {
        dictionary& bndLayersDict = meshDict_.subDict("boundaryLayers");
        if( bndLayersDict.isDict("patchBoundaryLayers") )
        {
            dictionary& patchBndLayers =
                bndLayersDict.subDict("patchBoundaryLayers");

            const wordList patchLayers = patchBndLayers.toc();

            forAll(patchLayers, patchI)
            {
                const word& pName = patchLayers[patchI];

                dictionary dict = patchBndLayers.subDict(pName);

                const std::map<word, wordList>::const_iterator it =
                    patchesFromPatch.find(pName);

                //- patch name may be a regex
                if( it != patchesFromPatch.end() )
                {
                    const wordList& newNames = it->second;

                    forAll(newNames, i)
                    {
                        patchBndLayers.add(newNames[i], dict);
                    }

                    patchBndLayers.remove(pName);
                }
            }
        }
    }
}

void checkMeshDict::updateRenameBoundary
(
    const std::map<word, wordList>& patchesFromPatch,
    const std::map<word, word>& patchTypes
)
{
    dictionary newDict;

    newDict.add("newPatchNames", dictionary());

    if( meshDict_.found("renameBoundary") )
    {
        const dictionary& dict = meshDict_.subDict("renameBoundary");

        //- transfer or generate the default name entry
        if( dict.found("defaultName") )
        {
            const word name(dict.lookup("defaultName"));
            newDict.add("defaultName", name);
        }

        //- transfer or generate the defaultType entry
        if( dict.found("defaultType") )
        {
            const word type(dict.lookup("defaultType"));
            newDict.add("defaultType", type);
        }

        if( dict.found("newPatchNames") )
        {
            //- stores the updated dictionary
            dictionary& newPatchesDict = newDict.subDict("newPatchNames");

            if( dict.isDict("newPatchNames") )
            {
                //- current state of the dictionary
                const dictionary& patchDicts = dict.subDict("newPatchNames");

                std::map<word, wordList>::const_iterator it;
                for(it=patchesFromPatch.begin();it!=patchesFromPatch.end();++it)
                {
                    const word& pName = it->first;
                    const wordList& newNames = it->second;

                    if( patchDicts.found(pName) )
                    {
                        //- patch renaming is already requested by the user
                        //- use the new name for all newly created patches
                        const dictionary& patchDict = patchDicts.subDict(pName);
                        if( !patchDict.found("newName") )
                            continue;
                        if( !patchDict.found("type") )
                            continue;

                        const word newName(patchDict.lookup("newName"));
                        const word newType(patchDict.lookup("type"));

                        forAll(newNames, i)
                        {
                            dictionary newPatchDict;
                            newPatchDict.add("newName", newName);
                            newPatchDict.add("type", newType);

                            newPatchesDict.add(newNames[i], newPatchDict);
                        }
                    }
                    else
                    {
                        //- rename all newly create patches
                        //- with the original name
                        forAll(newNames, i)
                        {
                            dictionary newPatchDict;

                            newPatchDict.add("newName", it->first);
                            std::map<word, word>::const_iterator tIter =
                                patchTypes.find(it->first);
                            newPatchDict.add("type", tIter->second);

                            newPatchesDict.add(newNames[i], newPatchDict);
                        }
                    }
                }
            }
            else
            {
                const PtrList<entry> patchEntries(dict.lookup("newPatchNames"));

                forAll(patchEntries, entryI)
                {
                    const word& pName = patchEntries[entryI].keyword();
                    dictionary patchDict(patchEntries[entryI].dict());

                    std::map<word, wordList>::const_iterator it =
                        patchesFromPatch.find(pName);

                    if( it == patchesFromPatch.end() )
                        continue;

                    const wordList& newNames = it->second;

                    forAll(newNames, i)
                        newPatchesDict.add(newNames[i], patchDict, true);
                }

                std::map<word, wordList>::const_iterator it;
                for(it=patchesFromPatch.begin();it!=patchesFromPatch.end();++it)
                {
                    const word& pName = it->first;
                    const wordList& newNames = it->second;

                    if( newPatchesDict.found(pName) )
                        continue;

                    //- rename all newly created patches
                    //- with the original name
                    forAll(newNames, i)
                    {
                        dictionary newPatchDict;

                        newPatchDict.add("newName", it->first);
                        std::map<word, word>::const_iterator tIter =
                            patchTypes.find(it->first);
                        newPatchDict.add("type", tIter->second);

                        newPatchesDict.add(newNames[i], newPatchDict);
                    }
                }
            }
        }
        else
        {
            //- newPatchNames is not used
            dictionary& newPatchesDict = newDict.subDict("newPatchNames");

            std::map<word, wordList>::const_iterator it;
            for(it=patchesFromPatch.begin();it!=patchesFromPatch.end();++it)
            {
                const wordList& newPatchNames = it->second;

                forAll(newPatchNames, i)
                {
                    const word& pName = newPatchNames[i];
                    dictionary newPatchDict;
                    newPatchDict.add("newName", it->first);
                    std::map<word, word>::const_iterator tIter =
                        patchTypes.find(it->first);
                    newPatchDict.add("type", tIter->second);

                    newPatchesDict.add(pName, newPatchDict);
                }
            }
        }

        //- delete all previus entries from the dictionary
        meshDict_.subDict("renameBoundary").clear();
    }
    else
    {
        //- create the dictionary if it has not existed before
        dictionary& newPatchesDict = newDict.subDict("newPatchNames");

        std::map<word, wordList>::const_iterator it;
        for(it=patchesFromPatch.begin();it!=patchesFromPatch.end();++it)
        {
            const wordList& newPatchNames = it->second;

            forAll(newPatchNames, i)
            {
                const word& pName = newPatchNames[i];
                dictionary newPatchDict;
                newPatchDict.add("newName", it->first);
                std::map<word, word>::const_iterator tIter =
                    patchTypes.find(it->first);
                newPatchDict.add("type", tIter->second);

                newPatchesDict.add(pName, newPatchDict);
            }
        }
    }

    meshDict_.add("renameBoundary", newDict, true);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

checkMeshDict::checkMeshDict
(
    IOdictionary& meshDict
)
:
    meshDict_(meshDict)
{
    checkEntries();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

checkMeshDict::~checkMeshDict()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkMeshDict::updateDictionaries
(
    const std::map<word, wordList>& patchesFromPatch,
    const std::map<word, word>& patchTypes,
    const bool renamePatches
)
{
    updatePatchCellSize(patchesFromPatch);

    updateSubsetCellSize(patchesFromPatch);

    updateLocalRefinement(patchesFromPatch);

    updateKeepCellsIntersectingPatches(patchesFromPatch);

    updateRemoveCellsIntersectingPatches(patchesFromPatch);

    updateObjectRefinements(patchesFromPatch);

    updateBoundaryLayers(patchesFromPatch);

    if( renamePatches )
        updateRenameBoundary(patchesFromPatch, patchTypes);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
