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
    refineBoundaryLayers

Description
    Refine existing boundary layers

SourceFiles
    refineBoundaryLayers.C
    refineBoundaryLayersFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef refineBoundaryLayers_H
#define refineBoundaryLayers_H

#include "objectRegistry.H"
#include "Time.H"
#include "polyMeshGenModifier.H"
#include "meshSurfaceEngine.H"
#include "DynList.H"
#include "labelLongList.H"
#include "labelPair.H"

#include <map>
#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class refineBoundaryLayers Declaration
\*---------------------------------------------------------------------------*/

class refineBoundaryLayers
{
        //- Reference to the mesh
        polyMeshGen& mesh_;

        //- pointer to mesh surface engine
        mutable meshSurfaceEngine* msePtr_;

        //- global number of boundary layers
        label globalNumLayers_;

        //- global thickness ratio
        scalar globalThicknessRatio_;

        //- global maximum thickness of the first layer
        scalar globalMaxThicknessFirstLayer_;

        //- number of boundary layers for user-selected patches
        std::map<word, label> numLayersForPatch_;

        //- local thickness ratio for selected patches
        std::map<word, scalar> thicknessRatioForPatch_;

        //- local maximum layer thickness for selected patches
        std::map<word, scalar> maxThicknessForPatch_;

        //- allow discontinuous layers for patch
        std::set<word> discontinuousLayersForPatch_;

        //- check whether the refinement is already executed
        bool done_;

        //- a flag whether a 2D mesh generation is active or not
        bool is2DMesh_;

        //- information about existing boundary layers at patches
        //- only available layers
        labelList layerAtPatch_;

        //- which patches are part of a single layer
        List<DynList<word> > patchesInLayer_;

        //- a containing the number of layers which shall be generated above
        //- a boundary face
        labelList nLayersAtBndFace_;

        //- a list of edges which shall be refined
        LongList<edge> splitEdges_;

        //- split edges at point
        VRWGraph splitEdgesAtPoint_;

        //- new vertices for on edges which shall be refined
        VRWGraph newVerticesForSplitEdge_;

        //- a graph containing information which new faces were generated
        //- from an existing face
        VRWGraph facesFromFace_;

        //- a graph containing faces after layer refinement
        VRWGraph newFaces_;

    // Private member functions
        //- Return reference to meshSurfaceEngine
        const meshSurfaceEngine& surfaceEngine() const;

        //- analyse layers to check their topology
        void analyseLayers();

        //- calculate addressing for a boundary cell
        void calculateAddressing
        (
            const label bfI,
            label& baseFace,
            DynList<edge, 48>& edges,
            DynList<DynList<label, 2>, 48>& edgeFaces,
            DynList<DynList<label, 10>, 24>& faceEdges
        ) const;

        //- find bnd layer hairs for a boundary face
        bool findHairsForFace(const label, DynList<edge>&) const;

        //- find edges which shall be split due to refinement
        bool findSplitEdges();

        //- generate new points on edges, faces and in cells
        void generateNewVertices();

        //- refine a given face and return the new faces
        //- generates new points at cross-split faces
        void refineFace
        (
            const face& f,
            const FixedList<label, 2>& nLayersInDirection,
            DynList<DynList<label, 4>, 128>& newFaces
        );

        //- generate a matrix of points generated by splitting a face
        //- and return them in the local i, j system of the face
        void sortFacePoints
        (
            const label faceI,
            DynList<DynList<label> >& facePoints,
            const label transpose = false
        ) const;

        //- generate a matrix of faces generated by splitting a face
        //- and return them in the local i, j, system of the face
        void sortFaceFaces
        (
            const label faceI,
            DynList<DynList<label> >& faceFaces,
            const label transpose = false
        ) const;

        //- map split edges onto a cell
        void generateNewFaces();

        //- generate new cells for a prism with one boundary face
        void generateNewCellsPrism
        (
            const label cellI,
            DynList<DynList<DynList<label, 8>, 10> >& cellsFromCell
        );

        //- a helper function which stores faces generated from
        //- an existing face into new cells
        void storeFacesIntoCells
        (
            const label faceI,
            const bool reverseOrientation,
            const label normalDirection,
            const bool maxCoordinate,
            const label nLayersI,
            const label nLayersJ,
            const label nLayersK,
            DynList<DynList<DynList<label, 4>, 6>, 256>& cellsFromCell
        ) const;

        //- generate new cells and add them to the mesh
        void generateNewCells();

    // Nested classes
        class refineEdgeHexCell
        {
            // Private data
                //- label of cell
                const label cellI_;

                //- number of cells in local direction i
                label nLayersI_;

                //- number of cells in locatiol direction j
                label nLayersJ_;

                //- container for new cells
                DynList<DynList<DynList<label, 4>, 6>, 256> cellsFromCell_;

                //- reference to the boundary layer class
                refineBoundaryLayers& bndLayers_;

                //- faces sorted into directions of a hex shape
                FixedList<label, 6> faceInDirection_;

                //- information about orientation of faces
                //- false means the orientation as expected
                //- true means wrong orientation
                FixedList<bool, 6> faceOrientation_;

                //- points on cross-split faces
                FixedList<DynList<DynList<label> >, 2> cellPoints_;

            // Private member functions
                //- populate faceInDirection_nad wrongFaceOrientation_
                void determineFacesInDirections();

                //- populate new cells with new faces generated from already
                //- existing faces
                void populateExistingFaces();

                //- generate new internal faces and tore them to new cells
                void generateMissingFaces();

        public:

            // Constructor
                //- construct from cell label and the refineBoundaryLayers
                refineEdgeHexCell(const label cellI, refineBoundaryLayers& ref);

            // Public member functions
                inline const DynList<DynList<DynList<label, 4>, 6>, 256>&
                newCells() const
                {
                    return cellsFromCell_;
                }
        };

        class refineCornerHexCell
        {
            // Private data
                //- label of cell
                const label cellI_;

                //- number of cells in local direction i
                label nLayersI_;

                //- number of cells in local direction j
                label nLayersJ_;

                //- number of cells in local direction k
                label nLayersK_;

                //- split edge in directions
                FixedList<label, 3> splitEdgeInDirection_;

                //- container for new cells
                DynList<DynList<DynList<label, 4>, 6>, 256> cellsFromCell_;

                //- reference to the boundary layer class
                refineBoundaryLayers& bndLayers_;

                //- faces sorted into directions of a hex shape
                FixedList<label, 6> faceInDirection_;

                //- information about orientation of faces
                //- false means the orientation as expected
                //- true means wrong orientation
                FixedList<bool, 6> faceOrientation_;

                //- points on cross-split faces
                FixedList<DynList<DynList<label> >, 6> facePoints_;

                //- points inside the cell
                DynList<DynList<DynList<label> > > cellPoints_;

            // Private member functions
                //- populate faceInDirection_nad wrongFaceOrientation_
                void determineFacesInDirections();

                //- populate new cells with new faces generated from already
                //- existing faces
                void populateExistingFaces();

                //- generate missing points inside the cell
                void generateNewPoints();

                //- generate new internal faces and tore them to new cells
                void generateMissingFaces();

        public:

            // Constructor
                //- construct from cell label and the refineBoundaryLayers
                refineCornerHexCell
                (
                    const label cellI,
                    refineBoundaryLayers& ref
                );

            // Public member functions
                inline const DynList<DynList<DynList<label, 4>, 6>, 256>&
                newCells() const
                {
                    return cellsFromCell_;
                }
        };

    // Private member functions

        //- Disallow bitwise copy construct
        refineBoundaryLayers(const refineBoundaryLayers&);

        //- Disallow bitwise assignment
        void operator=(const refineBoundaryLayers&);

public:

    // Constructors

        //- Construct from mesh reference
        refineBoundaryLayers(polyMeshGen& mesh);

    // Destructor
        ~refineBoundaryLayers();

    // Public member functions
        //- set no refinement flag
        void avoidRefinement();

        //- activate 2D layer refinement
        void activate2DMode();

        //- set the global number of boundary layers
        void setGlobalNumberOfLayers(const label nLayers);

        //- set the global thickness ratio (default is 1)
        void setGlobalThicknessRatio(const scalar thicknessRatio);

        //- set the maximum thickness of the first boundary layer
        void setGlobalMaxThicknessOfFirstLayer(const scalar maxThickness);

        //- set the number of layers for a patch
        //- the settings override the global settings
        void setNumberOfLayersForPatch
        (
            const word& patchName,
            const label nLayers
        );

        //- set the thickness ratio for a patch
        //- it overrides the global settings
        void setThicknessRatioForPatch
        (
            const word& patchName,
            const scalar thicknessRatio
        );

        //- set the maximum thickness of the first layer for a patch
        void setMaxThicknessOfFirstLayerForPatch
        (
            const word& patchName,
            const scalar maxThickness
        );

        //- set whether the settings for a given patch are valid for the
        //- patch only, or whether they extend over th whole sheet
        //- the selected patch belongs to
        //- the default behaviour is to apply the patch settings to the whole
        //- sheet
        void setInteruptForPatch(const word& patchName);

        //- performs refinement based on the given settings
        void refineLayers();

        //- read the settings from dictionary
        static void readSettings(const dictionary&, refineBoundaryLayers&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
