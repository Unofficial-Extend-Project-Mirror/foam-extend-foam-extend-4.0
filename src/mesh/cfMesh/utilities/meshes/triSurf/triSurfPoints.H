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
    triSurfPoints

Description
    Point coordinates information for triangulated surface

SourceFiles
    triSurfPoints.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfPoints_H
#define triSurfPoints_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshSubset.H"
#include "pointField.H"
#include "DynList.H"
#include "Map.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class triSurfPoints
{
protected:

    // Protected data
        //- list of vertices
        pointField points_;

        //- map of point subsets
        Map<meshSubset> pointSubsets_;

    // Disallow bitwise assignment
        void operator=(const triSurfPoints&);

        triSurfPoints(const triSurfPoints&);

public:

    // Constructors
        //- Null constructor
        triSurfPoints();

        //- Construct from vertices
        triSurfPoints(const pointField& points);

    // Destructor
        ~triSurfPoints();

    // Member functions
        //- return the number of points
        inline label nPoints() const;

        //- access to points
        inline const pointField& points() const;

        //- append a vertex to the end of the list
        inline void appendVertex(const point& p);

        //- point subsets
        label addPointSubset(const word&);
        void removePointSubset(const label);
        word pointSubsetName(const label) const;
        label pointSubsetIndex(const word&) const;
        inline void addPointToSubset(const label, const label);
        inline void removePointFromSubset(const label, const label);
        inline void pointInSubsets(const label, DynList<label>&) const;
        inline void pointSubsetIndices(DynList<label>&) const;
        template<class ListType>
        inline void pointsInSubset(const label, ListType&) const;
        template<class ListType>
        inline void updatePointSubsets(const ListType&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triSurfPointsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
