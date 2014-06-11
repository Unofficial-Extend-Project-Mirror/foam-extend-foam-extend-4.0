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
    triSurfFeatureEdges

Description
    User-selected list of feature edges used in the meshing process

SourceFiles
    triSurfFeatureEdges.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfFeatureEdges_H
#define triSurfFeatureEdges_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "edgeLongList.H"
#include "meshSubset.H"
#include "Map.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class triSurfFeatureEdges
{
protected:

    // Protected data
        //- list of feature edges
        edgeLongList featureEdges_;

        //- map of edge subsets
        Map<meshSubset> featureEdgeSubsets_;

    // Disallow bitwise assignment
        void operator=(const triSurfFeatureEdges&);

        triSurfFeatureEdges(const triSurfFeatureEdges&);

public:

    // Constructors
        //- Null constructor
        triSurfFeatureEdges();

        //- Construct from feature edges
        triSurfFeatureEdges(const edgeLongList& triangles);

    // Destructor
        ~triSurfFeatureEdges();

    // Member functions
        //- return the number of feature edges
        inline label nFeatureEdges() const;

        //- access to feature edges
        inline const edgeLongList& featureEdges() const;

        //- append an edge to the end of the list
        inline void appendFeatureEdge(const edge& e);

        //- point subsets
        label addEdgeSubset(const word&);
        void removeEdgeSubset(const label);
        word edgeSubsetName(const label) const;
        label edgeSubsetIndex(const word&) const;
        inline void addEdgeToSubset(const label, const label);
        inline void removeEdgeFromSubset(const label, const label);
        inline void edgeInSubsets(const label, DynList<label>&) const;
        inline void edgeSubsetIndices(DynList<label>&) const;
        template<class ListType>
        inline void edgesInSubset(const label, ListType&) const;
        template<class ListType>
        inline void updateEdgeSubsets(const ListType&);

    // Operators
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triSurfFeatureEdgesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
