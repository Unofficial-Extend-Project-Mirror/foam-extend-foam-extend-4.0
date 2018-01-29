/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    Foam::polyhedralRefinement

Description
    Isotropic refinement of polyhedral cells using the mesh modifier engine

    Each polyhedral cell is split by the following procedure:
    1. Adding points at the edge centre, face centre and cell centre,
    2. Adding cells n cells where n is the number of points of the cell,
    3. Splitting each face into multiple faces going from:
       existing corner point -> new edge centre point -> new face centre
       point -> other new edge centre point (sharing the same corner point)
    4. Adding internal faces going from:
       new edge centre point -> new face centre point -> new cell centre
       point -> other new face centre point (sharing the same edge)

SourceFiles
    polyhedralRefinement.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

Notes
    Generalisation of hexRef8 for polyhedral cells and refactorisation using
    polyMesh modifier engine.

\*---------------------------------------------------------------------------*/

#ifndef polyhedralRefinement_H
#define polyhedralRefinement_H

#include "polyMeshModifier.H"
#include "labelIOList.H"
#include "removeFaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class polyhedralRefinement Declaration
\*---------------------------------------------------------------------------*/

class polyhedralRefinement
:
    public polyMeshModifier
{
public:

    // Public enumeration for refinement status
    enum refinementStatus
    {
        UNREFINED = -1,
        UNCHANGED = 0,
        REFINED = 1
    };


private:

    // Private data

        //- Reference to polyMesh for easy access in helper functions
        const polyMesh& mesh_;


        // Refinement control and handling

            //- List of cells to refine in this time step
            mutable labelList cellsToRefine_;

            //- List of split point labels to unrefine in this time step
            mutable labelList splitPointsToUnrefine_;

            //- Cell refinement level
            mutable labelIOList cellLevel_;

            //- Point refinement level
            mutable labelIOList pointLevel_;

            //- Helper list for original (old) cells that will be refined or
            //  unrefined. The list is updated in setPolyhedralRefinement and
            //  setPolyhedralUnfereinement and is used in updateMesh to update
            //  the cellLevel on the new mesh.
            //  Values stored in the list are:
            //  a) UNREFINED = -1 = cell is unrefined
            //  b) UNCHANGED =  0 = cell is untouched
            //  c) REFINED   = +1 = cell is refined
            mutable labelList refinementLevelIndicator_;

            //- Typical edge length between unrefined points
            scalar level0EdgeLength_;

            //- Face remover engine
            mutable removeFaces faceRemover_;

            //- Maximum number of cells in the mesh. Note: not strictly enforced
            label maxCells_;

            //- Maximum number of refinement levels for a given cell
            label maxRefinementLevel_;

            //- Switch whether to use point based consistency on refinement
            Switch pointBasedConsistency_;

            //- Number of buffer layers for refinement
            label nBufferLayers_;

            //- Flag indicating whether we currently have topo changes
            //  Set to true in setRefinement member function
            //  Set to false in updateMesh member function
            //  Note: must be reduced. Even if a certain processor does not go
            //  through topo changes, it must call updateMesh to correctly
            //  update point levels of points on coupled patches
            mutable bool changingTopology_;


    // Private Member Functions

        // Helper functions

            //- Get least cell level such that the face has at least three
            //  points smaller than the level
            label getAnchorLevel(const label faceI) const;

            //- Calculate level0EdgeLength_ (constructor helper)
            void calcLevel0EdgeLength();

            //- Set file instance for cellLevel_ and pointLevel_
            void setInstance(const fileName& inst) const;

            //- Extend cells to refine. Given a markup field of refinement
            //  candidates (true for cells to refine), marks the neigbhours of
            //  marked cells as well
            void extendMarkedCells(boolList& refineCell) const;


        // Global topology modification functions (operate on whole polyMesh)
        // and directly use local topology modification functions below

            //- Refine polyhedral cells. All cellsToRefine_ cells will be split
            //  into n cells where n is the number of points per cell.
            //  Guarantees that the 0th element is the original cell label.
            //  Mapping:
            //  -split cells: n new ones get added from original
            //  -split faces: original gets modified; n - 1 new ones get added
            //                from original
            //  -added internal faces: added from original cell face (if
            //   that was internal) or created out-of-nothing (so will not
            //   get mapped!). Note: could make this inflate from point but
            //   that will allocate interpolation.
            //  -points added to split edge: added from edge start()
            //  -midpoints added: added from cellPoints[0].
            void setPolyhedralRefinement(polyTopoChange& ref) const;

            //- Remove some of the previously performed refinement. Uses
            //  splitPointLabels_ to determine the unrefinement.
            //  All n pointCells of a split point will be combined into
            //  the lowest numbered cell of those n.
            void setPolyhedralUnrefinement(polyTopoChange& ref) const;


        // Local topology modification functions (operate on cells/faces)

            //- Adds a face on top of existing faceI. Reverses if nessecary
            label addFace
            (
                polyTopoChange& ref,
                const label faceI,
                const face& newFace,
                const label own,
                const label nei
            ) const;

            //- Adds internal face from point. No checks on reversal
            label addInternalFace
            (
                polyTopoChange& ref,
                const label meshFaceI,
                const label meshPointI,
                const face& newFace,
                const label own,
                const label nei
            ) const;

            //- Modifies existing faceI for either new owner/neighbour or new face
            //  points. Reverses if nessecary
            void modifyFace
            (
                polyTopoChange& ref,
                const label faceI,
                const face& newFace,
                const label own,
                const label nei
            ) const;

            //- Create all internal faces of split cellI into n cells where n is the
            //  number of cell points
            void createInternalFaces
            (
                const labelListList& cellAnchorPoints,
                const labelListList& cellAddedCells,
                const labelList& cellMidPoint,
                const labelList& faceMidPoint,
                const labelList& faceAnchorLevel,
                const labelList& edgeMidPoint,
                const label cellI,
                polyTopoChange& ref
            ) const;


        // Topological change helper functions

            //- Get cell added to point of cellI (if any)
            label getAnchorCell
            (
                const labelListList& cellAnchorPoints,
                const labelListList& cellAddedCells,
                const label cellI,
                const label faceI,
                const label pointI
            ) const;

            //- Set new owner and neighbour (in unspecified order) of pointI
            //  on faceI
            void setNewFaceNeighbours
            (
                const labelListList& cellAnchorPoints,
                const labelListList& cellAddedCells,
                const label faceI,
                const label pointI,

                label& own,
                label& nei
            ) const;

            //- Store vertices from startFp up to face split point.
            //  Used when splitting face into n faces where n is the number of
            //  points in a face (or number of edges)
            void walkFaceToMid
            (
                const labelList& edgeMidPoint,
                const label cLevel,
                const label faceI,
                const label startFp,
                dynamicLabelList& faceVerts
            ) const;

            //- Same as walkFaceToMid but now walk back
            void walkFaceFromMid
            (
                const labelList& edgeMidPoint,
                const label cLevel,
                const label faceI,
                const label startFp,
                dynamicLabelList& faceVerts
            ) const;

            //- Get index of point with minimum point level
            label findMinLevel(const labelList& f) const;

            //- Get index of point with maximum point level
            label findMaxLevel(const labelList& f) const;

            //- Count number of vertices <= anchorLevel for a given face
            label countAnchors
            (
                const labelList& f,
                const label anchorLevel
            ) const;

            //- Find index of point with wantedLevel, starting from fp
            label findLevel
            (
                const face& f,
                const label startFp,
                const bool searchForward,
                const label wantedLevel
            ) const;

            //- Store in maps correspondence from midpoint to anchors and
            //  faces. Used when creating internal faces
            label storeMidPointInfo
            (
                const labelListList& cellAnchorPoints,
                const labelListList& cellAddedCells,
                const labelList& cellMidPoint,
                const labelList& edgeMidPoint,
                const label cellI,
                const label faceI,
                const bool faceOrder,
                const label midPointI,
                const label anchorPointI,
                const label faceMidPointI,

                Map<edge>& midPointToAnchors,
                Map<edge>& midPointToFaceMids,
                polyTopoChange& ref
            ) const;

            //- If p0 and p1 are existing vertices check if edge is split and insert
            //  splitPoint. Used with storing mid point
            void insertEdgeSplit
            (
                const labelList& edgeMidPoint,
                const label p0,
                const label p1,
                dynamicLabelList& verts
            ) const;


        // Debug functions

            //- Check orientation of added internal face
            void checkInternalOrientation
            (
                polyTopoChange& ref,
                const label cellI,
                const label faceI,
                const point& ownPt,
                const point& neiPt,
                const face& newFace
            ) const;

            //- Check orientation of a new boundary face
            void checkBoundaryOrientation
            (
                polyTopoChange& ref,
                const label cellI,
                const label faceI,
                const point& ownPt,
                const point& boundaryPt,
                const face& newFace
            ) const;


        // Refinement/unrefinement consistency checks

            //- Updates cellsToRefine such that a face consistent 2:1 refinement
            //  is obtained. Returns local number of cells changed
            label faceConsistentRefinement(boolList& cellsToRefine) const;

            //- Updates cellsToRefine such that a point consistent 4:1 refinement
            //  is obtained. Returns local number of cells changed
            label pointConsistentRefinement(boolList& cellsToRefine) const;

            //- Updates cellsToUnrefine such that a face consistent 2:1
            //  unrefinement is obtained. Returns local number of points changed
            label faceConsistentUnrefinement(boolList& cellsToUnrefine) const;

            //- Updates cellsToUnrefine such that a point consistent 4:1
            //  unrefinement is obtained. Returns local number of points changed
            label pointConsistentUnrefinement
            (
                const boolList& splitPointsToUnrefine,
                boolList& cellsToUnrefine
            ) const;


        // Copy control

            //- Disallow default bitwise copy construct
            polyhedralRefinement(const polyhedralRefinement&);

            //- Disallow default bitwise assignment
            void operator=(const polyhedralRefinement&);


public:

    //- Runtime type information
    TypeName("polyhedralRefinement");


    // Constructors

        //- Construct from dictionary
        polyhedralRefinement
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyTopoChanger& mme
        );


    // Destructor

        virtual ~polyhedralRefinement();


    // Member Functions

        // Access

            const labelIOList& cellLevel() const
            {
                return cellLevel_;
            }

            const labelIOList& pointLevel() const
            {
                return pointLevel_;
            }

            //- Typical edge length between unrefined points
            scalar level0EdgeLength() const
            {
                return level0EdgeLength_;
            }


        // Edit

            //- Set cells to refine given a list of refinement
            //  candidates. Refinement candidates are extended within the
            //  function due to possible 4:1 conflicts and specified number of
            //  buffer layers.
            //  Note: must be called BEFORE setSplitPointsToUnrefine
            void setCellsToRefine(const labelList& refinementCellCandidates);

            //- Set split points to unrefine given a list of all mesh points
            //  that are candidates for unrefinement. Split points are
            //  determined as a subset of unrefinement candidates, avoiding
            //  splitting points of cells that are going to be refined at the
            //  same time and ensuring consistent unrefinement.
            //  Note: must be called AFTER setCellsToRefine
            void setSplitPointsToUnrefine
            (
                const labelList& unrefinementPointCandidates
            );


        // Inherited interface from polyMeshModifier

            //- Check for topology change
            virtual bool changeTopology() const;

            //- Insert the polyhedral refinement/unrefinement into the
            //  topological change
            virtual void setRefinement(polyTopoChange&) const;

            //- Modify motion points to comply with the topological change
            virtual void modifyMotionPoints(pointField& motionPoints) const;

            //- Force recalculation of locally stored data on topological change
            virtual void updateMesh(const mapPolyMesh&);

            //- Write
            virtual void write(Ostream&) const;

            //- Write dictionary
            virtual void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
