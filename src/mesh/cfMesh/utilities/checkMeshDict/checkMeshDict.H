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

Class
    checkMeshDict

Description
    Check whether the meshDict file is set correctly

SourceFiles
    checkMeshDict.C

\*---------------------------------------------------------------------------*/

#ifndef checkMeshDict_H
#define checkMeshDict_H

#include "objectRegistry.H"
#include "Time.H"
#include "IOdictionary.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class checkMeshDict Declaration
\*---------------------------------------------------------------------------*/

class checkMeshDict
{
        //- Reference to the mesh
        IOdictionary& meshDict_;

    // Private member functions
        //- check settings for cell size in meshDict
        void checkBasicSettings() const;

        //- check patchCellSize entry
        void checkPatchCellSize() const;

        //- check subsetCellSize entry
        void checkSubsetCellSize() const;

        //- check local refinement level
        void checkLocalRefinementLevel() const;

        //- check keepCellsIntersectingPatches entry
        void checkKeepCellsIntersectingPatches() const;

        //- check removeCellsIntersectingPatches entry
        void checkRemoveCellsIntersectingPatches() const;

        //- check objectRefinements entry
        void checkObjectRefinements() const;

        //- check anisotropic sources
        void checkAnisotropicSources() const;

        //- check surfaceRefinements entry
        void checkSurfaceRefinements() const;

        //- check edgeMeshRefinements entry
        void checkEdgeMeshRefinements() const;

        //- check entry for boundary layers
        void checkBoundaryLayers() const;

        //- check renameBoundary entry
        void checkRenameBoundary() const;

        //- perform all checks
        void checkEntries() const;

        //- update patchCellSize entry
        void updatePatchCellSize(const std::map<word, wordList>&);

        //- update subsetCellSize entry
        void updateSubsetCellSize(const std::map<word, wordList>&);

        //- update local refinement
        void updateLocalRefinement(const std::map<word, wordList>&);

        //- check keepCellsIntersectingPatches entry
        void updateKeepCellsIntersectingPatches
        (
            const std::map<word, wordList>&
        );

        //- check removeCellsIntersectingPatches entry
        void updateRemoveCellsIntersectingPatches
        (
            const std::map<word, wordList>&
        );

        //- check objectRefinements entry
        void updateObjectRefinements(const std::map<word, wordList>&);

        //- check entry for boundary layers
        void updateBoundaryLayers(const std::map<word, wordList>&);

        //- check renameBoundary entry
        void updateRenameBoundary
        (
            const std::map<word, wordList>&,
            const std::map<word, word>&
        );

    public:

    // Constructors

        //- Construct from IOdictionary
        checkMeshDict(IOdictionary& meshDict);

    // Destructor
        ~checkMeshDict();

    // Public member functions

        //- update meshDict based on modification of patches in the surface
        void updateDictionaries
        (
            const std::map<word, wordList>& patchesForPatch,
            const std::map<word, word>& patchTypes,
            const bool renamePatches = true
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
