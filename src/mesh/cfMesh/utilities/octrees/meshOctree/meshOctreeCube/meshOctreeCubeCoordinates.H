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
    meshOctreeCubeCoordinates

Description
    A cube stores information needed for mesh generation

SourceFiles
    meshOctreeCubeCoordinates.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeCubeCoordinates_H
#define meshOctreeCubeCoordinates_H

#include "label.H"
#include "direction.H"
#include "FixedList.H"
#include "boundBox.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class triSurf;
class boundBox;

/*---------------------------------------------------------------------------*\
                Class meshOctreeCubeCoordinates Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeCubeCoordinates
{
    // Private data
        //- coordinates in the octree structure
        label posX_;
        label posY_;
        label posZ_;

        //- cube level in the octree structure
        direction level_;

public:

    // Static data
        //- edge nodes for an octree cube
        static const label edgeNodes_[12][2];

        //- cube nodes making each face
        static const label faceNodes_[6][4];

        //- node-faces addressing for the cube
        static const label nodeFaces_[8][3];

        //- face-edges addressing for the octree cube
        static const label faceEdges_[6][4];

        //- edge-faces addressing for the octree cube
        static const label edgeFaces_[12][2];

        //- return the opposite face of each cube face
        static const label oppositeFace_[6];

    // Constructors
        //- Null constructor
        inline meshOctreeCubeCoordinates();

        //- Construct from coordinates and level
        explicit inline meshOctreeCubeCoordinates
        (
            const label posX,
            const label posY,
            const label posZ,
            const direction level
        );

        //- copy constructor
        inline meshOctreeCubeCoordinates
        (
            const meshOctreeCubeCoordinates&
        );

    // Destructor

        ~meshOctreeCubeCoordinates()
        {}

    // Member functions
        //- return level
        inline direction level() const;

        //- return x, y, z coordinates
        inline label posX() const;
        inline label posY() const;
        inline label posZ() const;

        //- return the coordinates of child cube at the given position
        inline meshOctreeCubeCoordinates refineForPosition(const label) const;

        //- return the coordinates of the parent at the level reduced by
        //- the given number of levels
        inline meshOctreeCubeCoordinates reduceLevelBy
        (
            const direction diff
        ) const;

        //- return the coordinates of the parent at the given level
        inline meshOctreeCubeCoordinates reduceToLevel(const direction) const;

        //- return the minimal coordinates of the child at the given level
        inline meshOctreeCubeCoordinates increaseToLevelMin
        (
            const direction l
        ) const;

        //- return the maximal coordinates of the child at the given level
        inline meshOctreeCubeCoordinates increaseToLevelMax
        (
            const direction l
        ) const;

        //- return min and max points
        inline void cubeBox(const boundBox&, point&, point&) const;

        //- calculate vertices
        void vertices(const boundBox&, FixedList<point, 8>&) const;

        //- return centre
        inline point centre(const boundBox&) const;

        //- return size
        inline scalar size(const boundBox&) const;

        //- edges of the cube
        void edgeVertices
        (
            const boundBox&,
            FixedList<FixedList<point, 2>, 12>&
        ) const;

        //- check if the surface triangle intersects the cube
        bool intersectsTriangle
        (
            const triSurf&,
            const boundBox&,
            const label
        ) const;

        bool intersectsTriangleExact
        (
            const triSurf&,
            const boundBox&,
            const label
        ) const;

        //- is a vertex inside the cube
        bool isVertexInside(const boundBox&, const point&) const;
        bool isPositionInside(const meshOctreeCubeCoordinates&) const;

        //- return the possible range of neighbour boxes one level
        //- higher than the current cube. It is useful for parallel octree
        inline void neighbourRange
        (
            meshOctreeCubeCoordinates& minCoord,
            meshOctreeCubeCoordinates& maxCoord
        ) const;

        //- check if the cube intersects a line
        bool intersectsLine(const boundBox&, const point&, const point&) const;

    // Operators

        inline meshOctreeCubeCoordinates operator+
        (
            const meshOctreeCubeCoordinates&
        ) const;
        inline void operator=(const meshOctreeCubeCoordinates&);
        inline bool operator==(const meshOctreeCubeCoordinates&) const;
        inline bool operator!=(const meshOctreeCubeCoordinates&) const;

        //- comparison of Morton codes
        inline bool operator<=(const meshOctreeCubeCoordinates&) const;
        inline bool operator>=(const meshOctreeCubeCoordinates&) const;

        friend inline Istream& operator>>
        (
            Istream&,
            meshOctreeCubeCoordinates&
        );
        friend inline Ostream& operator<<
        (
            Ostream&,
            const meshOctreeCubeCoordinates&
        );
};

//- Specify data associated with meshOctreeCubeCoordinates type is contiguous
template<>
inline bool contiguous<meshOctreeCubeCoordinates>() {return true;}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshOctreeCubeCoordinatesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
