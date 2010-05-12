/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::treeBoundBox

Description
    Standard boundBox + extra functionality for use in octree.

    Numbering of corner points is according to octant numbering.

    On the back plane (z=0):

    @verbatim
        Y
        ^
        |
        +--------+
        |2      3|
        |        |
        |        |
        |        |
        |0      1|
        +--------+->X
    @endverbatim

    For the front plane add 4 to the point labels.


SourceFiles
    treeBoundBoxI.H
    treeBoundBox.C

\*---------------------------------------------------------------------------*/

#ifndef treeBoundBox_H
#define treeBoundBox_H

#include "boundBox.H"
#include "direction.H"
#include "pointField.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Random;

/*---------------------------------------------------------------------------*\
                           Class treeBoundBox Declaration
\*---------------------------------------------------------------------------*/

class treeBoundBox
:
    public boundBox
{

private:

        //- To initialise faceNormals.
        static FixedList<vector, 6> calcFaceNormals();

public:

    // Static data members

        static const treeBoundBox greatBox;

        //- Bits used for octant/point coding. Every octant/corner point
        //  is the combination of three faces.
        enum octantBit
        {
            RIGHTHALF   = 0x1 << 0,
            TOPHALF     = 0x1 << 1,
            FRONTHALF   = 0x1 << 2
        };

        //- Face codes
        enum faceId
        {
            LEFT   = 0,
            RIGHT  = 1,
            BOTTOM = 2,
            TOP    = 3,
            BACK   = 4,
            FRONT  = 5
        };

        //- Bits used for face coding
        enum faceBit
        {
            NOFACE    = 0,
            LEFTBIT   = 0x1 << LEFT,
            RIGHTBIT  = 0x1 << RIGHT,
            BOTTOMBIT = 0x1 << BOTTOM,
            TOPBIT    = 0x1 << TOP,
            BACKBIT   = 0x1 << BACK,
            FRONTBIT  = 0x1 << FRONT,
        };

        //- Edges codes.
        //  E01 = edge between 0 and 1.
        enum edgeId
        {
            E01 = 0,
            E13 = 1,
            E23 = 2,
            E02 = 3,

            E45 = 4,
            E57 = 5,
            E67 = 6,
            E46 = 7,

            E04 = 8,
            E15 = 9,
            E37 = 10,
            E26 = 11
        };

        //- Face to point addressing
        static const faceList faces;

        //- Edge to point addressing
        static const edgeList edges;

        //- Per face the unit normal
        static const FixedList<vector, 6> faceNormals;


        //- Face on which neighbour is
        static direction neighbourFaceBits(const label&);


    // Constructors

        //- Construct null setting points to zero
        inline treeBoundBox();

        //- Construct from components
        inline treeBoundBox(const point& min, const point& max);

        //- Construct from components
        inline treeBoundBox(const boundBox& bb);

        //- Construct as the bounding box of the given pointField. Local
        //  processor domain only (no reduce as in boundBox)
        treeBoundBox(const UList<point>& points);

        //- Construct as subset of points
        treeBoundBox(const UList<point>&, const labelList& meshPoints);

        //- Construct from Istream
        treeBoundBox(Istream&);


    // Member functions

        // Access

            //- Smallest of length,height,width
            inline scalar minDim() const;

            //- Largest of length,height,width
            inline scalar maxDim() const;

            //- Average of length,height,width
            inline scalar avgDim() const;

            //- Typical dimension length,height,width
            inline scalar typDim() const;

            //- vertex coordinates. In octant coding.
            pointField points() const;


        // Check

            //- Corner point given octant
            inline point corner(const direction) const;

            //- Calculates midpoint
            inline point mid() const;

            //- Sub box given by octant number. Midpoint calculated.
            treeBoundBox subBbox(const direction) const;

            //- Sub box given by octant number. Midpoint provided.
            treeBoundBox subBbox(const point& mid, const direction) const;

            //- Returns octant number given point. Midpoint calculated.
            inline direction subOctant
            (
                const point& sample
            ) const;

            //- Returns octant number given point. Midpoint provided.
            static inline direction subOctant
            (
                const point& mid,
                const point& sample
            );

            //- Returns octant number given point. Midpoint calculated.
            //  onEdge set if sample on edge of subOctant
            inline direction subOctant
            (
                const point& mid,
                bool& onEdge
            ) const;

            //- Returns octant number given point. Midpoint provided.
            //  onEdge set if sample on edge of subOctant
            static inline direction subOctant
            (
                const point& mid,
                const point& sample,
                bool& onEdge
            );

            //- Returns octant number given intersection. Midpoint provided.
            //  onEdge set if sample on edge of subOctant. If onEdge
            // the direction vector determines which octant to use
            // (acc. to which octant the sample would be if it were moved
            //  along dir)
            static inline direction subOctant
            (
                const point& mid,
                const vector& dir,
                const point& sample,
                bool& onEdge
            );

            //- Calculates optimal order to look for nearest to sample. First
            //  will be the octant containing the sample, second the octant
            //  with boundary nearest to the sample etc.
            inline void searchOrder
            (
	         const point& sample,
                 FixedList<direction, 8>& octantOrder
            ) const;

            //- Overlaps other boundingbox?
            inline bool overlaps(const treeBoundBox&) const;

            //- Overlaps boundingSphere (centre + sqr(radius))?
            bool overlaps(const point&, const scalar radiusSqr) const;

            //- Intersects other boundingbox?
            inline bool intersects(const treeBoundBox&) const;

            //- Intersects segment; set point to intersection position,
            //  return true if intersection found.
            // (intPt argument used during calculation even if not intersecting)
            bool intersects(const point&, const point&, point& intPt) const;

            //- fully contains bb
            bool contains(const treeBoundBox& bb) const;

            //- Contains point? (inside or on edge)
            inline bool contains(const point&) const;

            //- Contains point? (only inside)
            bool containsNarrow(const point&) const;

            //- Contains point inside or
            //  on edge and moving in direction dir would cause it to go
            //  inside.
            bool contains(const vector& dir, const point&) const;

            //- Position of point relative to bb
            direction posBits(const point&) const;

            //- Calculate nearest and furthest (to sample) vertex coords of
            //  bounding box
            void calcExtremities
            (
                const point& sample,
                point& nearest,
                point& furthest
            ) const;

            //- Returns distance sample to furthest away corner.
            scalar maxDist(const point& sample) const;

            //- Compare distance to point with other bounding box
            //  return:
            //  -1 : all vertices of my bounding box are nearer than any of
            //       other
            //  +1 : all vertices of my bounding box are further away than
            //       any of other
            //   0 : none of the above.
            label distanceCmp(const point&, const treeBoundBox& other) const;

            //- Return slightly wider bounding box
            //  Extends all dimensions with s*span*Random::scalar01()
            inline treeBoundBox extend(Random&, const scalar s) const;

    // Friend Operators

        friend bool operator==(const treeBoundBox&, const treeBoundBox&);
        friend bool operator!=(const treeBoundBox&, const treeBoundBox&);


    // IOstream operator

        friend Istream& operator>>(Istream&, treeBoundBox&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "treeBoundBoxI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
