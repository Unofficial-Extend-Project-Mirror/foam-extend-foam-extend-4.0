/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::removeFaces

Description
    Given list of faces to remove insert all the topology changes. Contains
    helper function to get consistent set of faces to remove.

SourceFiles
    removeFaces.C

Notes
    All classes (hexRef8, undoableMeshCutter) that have a private member object
    of this class (removeFaces) should be rewritten such that they are
    polyMeshModifier classes. In order to allow for smooth transition between
    using the directTopoChange engine towards polyTopoChange engine, the
    member functions that actually do the refinement in this class are templated
    on a Type which can either be polyTopoChange engine or directTopoChange
    since they provide the same interface. VV, 16/Jan/2018.

\*---------------------------------------------------------------------------*/

#ifndef removeFaces_H
#define removeFaces_H

#include "Pstream.H"
#include "HashSet.H"
#include "Map.H"
#include "boolList.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class directTopoChange;
class face;
class mapPolyMesh;
class mapDistributePolyMesh;

/*---------------------------------------------------------------------------*\
                           Class removeFaces Declaration
\*---------------------------------------------------------------------------*/

class removeFaces
{
    // Private data

        //- Const reference to mesh
        const polyMesh& mesh_;

        //- Cosine of angles between boundary faces. Boundary faces can be
        //  merged only if angle between faces is greater than minCos
        const scalar minCos_;


    // Private Member Functions

        //- Change elements in cellRegion that are oldRegion to newRegion.
        //  Recurses to cell neighbours.
        void changeCellRegion
        (
            const label cellI,
            const label oldRegion,
            const label newRegion,

            labelList& cellRegion
        ) const;

        //- Changes region of connected set of faces. Returns number of changed
        //  faces
        label changeFaceRegion
        (
            const labelList& cellRegion,
            const boolList& removedFace,
            const labelList& nFacesPerEdge,
            const label faceI,
            const label newRegion,

            labelList& faceRegion
        ) const;

        //- Find and return all affected faces (including faces marked for removal)
        Xfer<boolList> affectedFaces
        (
            const labelList& cellRegion,
            const labelList& cellRegionMaster,
            const labelList& facesToRemove,
            const labelHashSet& edgesToRemove,
            const labelHashSet& pointsToRemove
        ) const;


        // Topological changes

            //- Debug: write set of faces to file in obj format.
            static void writeOBJ
            (
                const indirectPrimitivePatch&,
                const fileName&
            );

            //- Merge faceLabels into single face
            template<class TopoChangeEngine>
            void mergeFaces
            (
                const labelList& cellRegion,
                const labelList& cellRegionMaster,
                const labelHashSet& pointsToRemove,
                const labelList& faceLabels,

                TopoChangeEngine& ref
            ) const;

            //- Return face with all pointsToRemove removed
            face filterFace
            (
                const labelHashSet& pointsToRemove,
                const label
            ) const;

            //- Wrapper for meshMod.modifyFace. Reverses face if own>nei.
            template<class TopoChangeEngine>
            void modifyFace
            (
                const face& f,
                const label masterFaceID,
                const label own,
                const label nei,
                const bool flipFaceFlux,
                const label newPatchID,
                const bool removeFromZone,
                const label zoneID,
                const bool zoneFlip,

                TopoChangeEngine& ref
            ) const;



        //- Disallow default bitwise copy construct
        removeFaces(const removeFaces&);

        //- Disallow default bitwise assignment
        void operator=(const removeFaces&);


public:

    //- Runtime type information
    ClassName("removeFaces");


    // Constructors

        //- Construct from mesh and min cos of angle for boundary faces
        //  to be considered aligned. Set to >= 1 to disable checking
        //  and always merge (if on same patch)
        removeFaces(const polyMesh&, const scalar minCos);


    // Member Functions

        //- Given a set of faces to pierce calculates:
        //  - region for connected cells
        //  - mastercell for each region. This is the lowest numbered cell
        //    of all cells that get merged.
        //  - new set of faces which contains input set + additional ones
        //    where cells on both sides would have same mastercell.
        //  Returns number of regions.
        label compatibleRemoves
        (
            const labelList& inPiercedFaces,

            labelList& cellRegion,
            labelList& cellRegionMaster,
            labelList& outPiercedFaces
        ) const;


        //- Play commands into TopoChangeEngine to remove faces.
        template<class TopoChangeEngine>
        void setRefinement
        (
            const labelList& piercedFaces,
            const labelList& cellRegion,
            const labelList& pointRegionMaster,
            labelList& cellRegionMaster,

            TopoChangeEngine& ref
        ) const;

        //- Force recalculation of locally stored data on topological
        //  change. Does nothing
        void updateMesh(const mapPolyMesh&)
        {}

        //- Force recalculation of locally stored data for mesh
        //  distribution. Does nothing.
        void distribute(const mapDistributePolyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "removeFacesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
