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
    triSurfAddressing

Description
    Provides topology addressing of the triangulated surface

SourceFiles
    triSurfAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfAddressing_H
#define triSurfAddressing_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VRWGraph.H"
#include "edgeLongList.H"
#include "triSurfFacets.H"
#include "labelledTri.H"
#include "DynList.H"

#include <map>

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class triSurfAddressing
{
    // Private data
        //- const reference to the points
        const pointField& points_;

    // Topology addressing
        //- const reference to the facets
        const LongList<labelledTri>& facets_;

        //- facets connected to a point
        mutable VRWGraph* pointFacetsPtr_;

        //- edges in the triangulation
        mutable edgeLongList* edgesPtr_;

        //- labels of edges in the triangles
        mutable VRWGraph* facetEdgesPtr_;

        //- labels of triangles connected to an edge
        mutable VRWGraph* edgeFacetsPtr_;

        //- labels of edges connected to a point
        mutable VRWGraph* pointEdgesPtr_;

        //- facets connected tp a facet via edges
        mutable VRWGraph* facetFacetsEdgesPtr_;

    // Geometry data
        //- point normals
        mutable vectorField* pointNormalsPtr_;

        //- face normals
        mutable vectorField* facetNormalsPtr_;

        //- face centres
        mutable vectorField* facetCentresPtr_;

    // Private member functions
        //- calculate point-facets addressing
        void calculatePointFacets() const;

        //- calculate edges, facet-edges and edge-facets addressing
        void calculateEdges() const;

        //- calculate facet-edges addresing
        void calculateFacetEdges() const;

        //- calculate edge-facets addressing
        void calculateEdgeFacets() const;

        //- calculate point-edges addressing
        void calculatePointEdges() const;

        //- calculate facet-faceys addressing
        void calculateFacetFacetsEdges() const;

        //- calculate point normals
        void calculatePointNormals() const;

        //- calculate normals of facets
        void calculateFacetNormals() const;

        //- calculate centres of facets
        void calculateFacetCentres() const;

    // Disallow bitwise assignment
        void operator=(const triSurfAddressing&);

        triSurfAddressing(const triSurfAddressing&);

public:

    // Constructors
        //- Construct from surface triangles
        triSurfAddressing
        (
            const pointField& points,
            const LongList<labelledTri>& triangles
        );

    // Destructor
        ~triSurfAddressing();

    // Member functions
        //- return point-facets addressing
        inline const VRWGraph& pointFacets() const;

        //- return edges
        inline const LongList<edge>& edges() const;

        //- return facet-edges addressing
        inline const VRWGraph& facetEdges() const;

        //- return edge-facets addressing
        inline const VRWGraph& edgeFacets() const;

        //- return point-edges addressing
        inline const VRWGraph& pointEdges() const;

        //- return facet-facets addressing
        inline const VRWGraph& facetFacets() const;

        //- return point normals
        inline const pointField& pointNormals() const;

        //- return normals of facets
        inline const vectorField& facetNormals() const;

        //- return centres of facets
        inline const vectorField& facetCentres() const;

        //- delete all data
        void clearAddressing();

        //- delete geometry data
        void clearGeometry();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triSurfAddressingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
