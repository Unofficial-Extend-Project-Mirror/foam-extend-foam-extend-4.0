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
    meshOps

Description
    Various utility functions that perform mesh-related operations.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    meshOpsI.H
    meshOps.C

\*---------------------------------------------------------------------------*/

#ifndef meshOps_H
#define meshOps_H

#include "Map.H"
#include "point.H"
#include "label.H"
#include "scalar.H"
#include "HashSet.H"
#include "cellList.H"
#include "edgeList.H"
#include "faceList.H"
#include "triPointRef.H"
#include "tetPointRef.H"
#include "vectorField.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                        Namespace meshOps Declaration
\*---------------------------------------------------------------------------*/

namespace meshOps
{
    // Method to find the common-edge between two faces.
    inline bool findCommonEdge
    (
        const label first,
        const label second,
        const UList<labelList>& faceEdges,
        label& common
    );

    // Method to find the interior/boundary faces
    // for an input quad-face and adjacent triangle-prism cell.
    inline void findPrismFaces
    (
        const label fIndex,
        const label cIndex,
        const UList<face>& faces,
        const UList<cell>& cells,
        const UList<label>& neighbour,
        FixedList<face,2>& bdyf,
        FixedList<label,2>& bidx,
        FixedList<face,2>& intf,
        FixedList<label,2>& iidx
    );

    // Utility method to build a hull of cells
    // connected to the edge (for 2D simplical meshes)
    inline void constructPrismHull
    (
        const label eIndex,
        const UList<face>& faces,
        const UList<cell>& cells,
        const UList<label>& owner,
        const UList<label>& neighbour,
        const UList<labelList>& edgeFaces,
        labelList& hullTriFaces,
        labelList& hullCells
    );

    // Utility method to build a hull of cells (and faces)
    // around an edge (for 3D simplical meshes)
    inline void constructHull
    (
        const label eIndex,
        const UList<face>& faces,
        const UList<edge>& edges,
        const UList<cell>& cells,
        const UList<label>& owner,
        const UList<label>& neighbour,
        const UList<labelList>& faceEdges,
        const UList<labelList>& edgeFaces,
        const labelList& hullVertices,
        labelList& hullEdges,
        labelList& hullFaces,
        labelList& hullCells,
        labelListList& ringEntities
    );

    // Utility method to find the isolated point
    // given two triangular faces.
    inline label findIsolatedPoint
    (
        const face& baseFace,
        const face& checkFace
    );

    // Method to find the isolated point on a triangular face
    // that doesn't lie on the specified edge
    inline void findIsolatedPoint
    (
        const face& f,
        const edge& e,
        label& ptIndex,
        label& nextPtIndex
    );

    // Given a face and cell index, find the apex point for a tet cell.
    inline label tetApexPoint
    (
        const label cIndex,
        const label fIndex,
        const UList<face>& faces,
        const UList<cell>& cells
    );

    // Given a cell index, find the centroid / volume of a cell.
    inline void cellCentreAndVolume
    (
        const label cIndex,
        const vectorField& points,
        const UList<face>& faces,
        const UList<cell>& cells,
        const UList<label>& owner,
        vector& centre,
        scalar& volume
    );

    // Determine whether a segment intersects a triangular face
    inline bool segmentTriFaceIntersection
    (
        const triPointRef& faceToCheck,
        const linePointRef& edgeToCheck,
        vector& intPoint
    );

    // Determine whether the particular point lies
    // inside the given triangular face
    inline bool pointInTriFace
    (
        const triPointRef& faceToCheck,
        const vector& checkPoint,
        bool testCoplanarity
    );

    // Dijkstra's algorithm for the shortest path problem
    inline bool Dijkstra
    (
        const Map<point>& points,
        const Map<edge>& edges,
        const label startPoint,
        const label endPoint,
        Map<label>& pi
    );

    // Method to insert a label between two labels in a list
    inline void insertLabel
    (
        const label newLabel,
        const label labelA,
        const label labelB,
        labelList& list
    );

    // Utility method to replace a label in a given list
    inline void replaceLabel
    (
         const label original,
         const label replacement,
         labelList& list
    );

    // Utility method to size-up the list to include an item
    template <class Type>
    inline void sizeUpList
    (
        const Type& item,
        List<Type>& list
    );

    // Utility method to size-down the list to remove an item
    template <class Type>
    inline void sizeDownList
    (
        const Type& item,
        List<Type>& list
    );

    // Remove an item at a particular index in the list
    template <class Type>
    inline void removeIndex
    (
        const label index,
        List<Type>& list
    );

    // Parallel send
    inline void pWrite
    (
        const label toID,
        const label& data,
        bool blocking = true
    );

    // Parallel send (for fixed lists)
    template <class Type, unsigned Size>
    inline void pWrite
    (
        const label toID,
        const FixedList<Type,Size>& data
    );

    // Parallel send (for lists)
    template <class Type>
    inline void pWrite
    (
        const label toID,
        const UList<Type>& data
    );

    // Parallel receive
    inline void pRead
    (
        const label fromID,
        label& data,
        bool blocking = true
    );

    // Parallel receive (for fixed lists)
    template <class Type, unsigned Size>
    inline void pRead
    (
        const label fromID,
        FixedList<Type,Size>& data
    );

    // Parallel receive (for lists)
    template <class Type>
    inline void pRead
    (
        const label fromID,
        UList<Type>& data
    );

    // Wait for buffer transfer completion.
    inline void waitForBuffers();

    // Select a list of elements from connectivity,
    // and output to a VTK format
    inline void writeVTK
    (
        const polyMesh& mesh,
        const word& name,
        const labelList& cList,
        const label primitiveType,
        const UList<point>& meshPoints,
        const UList<edge>& edges = List<edge>(0),
        const UList<face>& faces = List<face>(0),
        const UList<cell>& cells = List<cell>(0),
        const UList<label>& owner = List<label>(0),
        const UList<scalar>& scalField = UList<scalar>(),
        const UList<label>& lablField = UList<label>(),
        const UList<vector>& vectField = UList<vector>()
    );

    // Actual routine to write out the VTK file
    inline void writeVTK
    (
        const polyMesh& mesh,
        const word& name,
        const label nPoints,
        const label nCells,
        const label nTotalCells,
        const vectorField& points,
        const labelListList& cpList = labelListList(0),
        const label primitiveType = 3,
        const Map<label>& reversePointMap = Map<label>(),
        const Map<label>& reverseCellMap = Map<label>(),
        const UList<scalar>& scalField = UList<scalar>(),
        const UList<label>& lablField = UList<label>(),
        const UList<vector>& vectField = UList<vector>()
    );

} // End namespace meshOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add additional operations missing in ops.H
template<class T>
class lessOp
{ public: T operator()(const T& x, const T& y) const { return x < y; } };

template<class T>
class lessEqOp
{ public: T operator()(const T& x, const T& y) const { return x <= y; } };

template<class T>
class greaterOp
{ public: T operator()(const T& x, const T& y) const { return x > y; } };

template<class T>
class greaterEqOp
{ public: T operator()(const T& x, const T& y) const { return x >= y; } };

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshOpsI.H"

#ifdef NoRepository
#    include "meshOps.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
