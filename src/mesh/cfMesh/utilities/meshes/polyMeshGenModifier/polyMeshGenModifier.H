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
    polyMeshGenModifier

Description
    Modifier for polyMeshGen

SourceFiles
    polyMeshGenModifier.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenModifier_H
#define polyMeshGenModifier_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGen.H"
#include "boolList.H"
#include "VRWGraph.H"
#include "demandDrivenData.H"

namespace Foam
{

// Forward declarations
class VRWGraphList;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGenModifier
{
    // Private data
        //- reference to the mesh
        polyMeshGen& mesh_;

        //- helper data for adding cells
        VRWGraph* pointFacesPtr_;

    // Private member functions
        //- calculate and return point faces
        inline void calculatePointFaces()
        {
            const faceListPMG& faces = mesh_.faces();

            pointFacesPtr_ = new VRWGraph();
            VRWGraph& pointFaces = *pointFacesPtr_;

            pointFaces.reverseAddressing(mesh_.points().size(), faces);
        };

        //- re-order positions of processor boundary faces
        //- they should comea immediately after the internal faces
        void reorderProcBoundaryFaces();

protected:

        VRWGraph& pointFaces()
        {
            if( !pointFacesPtr_ )
                calculatePointFaces();

            return *pointFacesPtr_;
        };

public:

    // Constructors
        //- Construct from the reference to the mesh
        polyMeshGenModifier(polyMeshGen& mesh)
        :
            mesh_(mesh),
            pointFacesPtr_(NULL)
        {
            //mesh_.clearOut();
        };

    // Destructor
        ~polyMeshGenModifier()
        {
            this->clearOut();
        };

    // Member functions
        //- access to mesh points
        inline pointFieldPMG& pointsAccess()
        {
            return mesh_.points_;
        };

        //- access to mesh faces
        inline faceListPMG& facesAccess()
        {
            return mesh_.faces_;
        };

        //- access to cells
        inline cellListPMG& cellsAccess()
        {
            return mesh_.cells_;
        };

        //- access to processor boundary data
        inline PtrList<processorBoundaryPatch>& procBoundariesAccess()
        {
            return mesh_.procBoundaries_;
        }

        //- access to boundary data
        inline PtrList<boundaryPatch>& boundariesAccess()
        {
            return mesh_.boundaries_;
        }

        //- functions which change the mesh
        //- reorder boundary faces
        void reorderBoundaryFaces();

        //- remove unused vertices
        void removeUnusedVertices();

        //- remove faces
        void removeFaces(const boolList& removeFace);

        //- remove duplicate faces from the mesh
        void removeDuplicateFaces();

        //- remove cells
        void removeCells
        (
            const boolList& removeCell,
            const bool removeProcFaces = true
        );

        //- add cells (vertices must be added)
        void addCells(const LongList<faceList>& cellFaces);
        void addCells(const VRWGraphList& cellFaces);
        void addCell(const faceList& cellFaces);

        //- replace the boundary with new boundary faces
        void replaceBoundary
        (
            const wordList& patchNames,
            const VRWGraph& boundaryFaces,
            const labelLongList& faceOwners,
            const labelLongList& facePatches
        );

        //- add additional faces into processor patches
        void addProcessorFaces
        (
            const VRWGraph& procFaces,
            const labelLongList& facePatches
        );

        //- add new processor patch and return its label
        label addProcessorPatch(const label otherProcLabel);

        //- remove empty processor patch
        bool removeEmptyProcessorPatches();

        //- add buffer cells needed for exporting the mesh in the format
        //- required by some solvers
        void addBufferCells();

        //- zip up topologically open cells
        void zipUpCells();

        //- reorder the cells and faces to reduce the matrix bandwidth
        void renumberMesh();

        //- clear out unnecessary data (pointFacesPtr_);
        inline void clearOut()
        {
            deleteDemandDrivenData(pointFacesPtr_);
        }

        //- clear out all allocated data
        inline void clearAll()
        {
            clearOut();
            mesh_.clearOut();
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
