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
    triSurfaceClassifyEdges

Description
    Divides the surface mesh into regions bounded by feature edges

SourceFiles
    triSurfaceClassifyEdges.C
    triSurfaceClassifyEdgesFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceClassifyEdges_H
#define triSurfaceClassifyEdges_H

#include "triSurf.H"
#include "VRWGraph.H"
#include "boolList.H"
#include "direction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class meshOctree;

/*---------------------------------------------------------------------------*\
                Class triSurfaceClassifyEdges Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceClassifyEdges
{
    // Private data
        //- reference to meshOctree
        const meshOctree& octree_;

        //- flags for surface edges
        List<direction> edgeTypes_;

        //- orientation of facet's normal (0 - unknown, 1- outward, 2- inward)
        List<direction> facetOrientation_;

    // Private member functions
        //- check the orientation of the patches in the triangulated surface
        void checkOrientation();

        //- classify edges based on the orientation of the surface facets
        void classifyEdgesTypes();

        //- Disallow default bitwise copy construct
        triSurfaceClassifyEdges(const triSurfaceClassifyEdges&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceClassifyEdges&);

public:

    // Enumerators

        enum edgeType_
        {
            NONE = 0,
            CONVEXEDGE = 1,
            CONCAVEEDGE = 2,
            FLATSURFACEEDGE = 4,
            FEATUREEDGE = 8
        };

    // Constructors

        //- Construct from octree
        triSurfaceClassifyEdges(const meshOctree& octree);

    // Destructor

        ~triSurfaceClassifyEdges();


    // Member Functions

        //- return the edge type according to the above enumerator
        inline direction edgeType(const label edgeI) const
        {
            return edgeTypes_[edgeI];
        }

        //- return the list of edge classification
        const List<direction>& edgeTypes() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
