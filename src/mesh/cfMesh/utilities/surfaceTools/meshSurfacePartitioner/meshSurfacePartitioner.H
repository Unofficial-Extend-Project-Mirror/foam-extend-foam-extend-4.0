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
    meshSurfacePartitioner

Description
    Finds corners and edge points at the surface of the volume mesh

SourceFiles
    meshSurfacePartitioner.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfacePartitioner_H
#define meshSurfacePartitioner_H

#include "meshSurfaceEngine.H"
#include "VRWGraph.H"
#include "HashSet.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class meshSurfacePartitioner Declaration
\*---------------------------------------------------------------------------*/

class meshSurfacePartitioner
{
    // Private data
        //- reference to mesh surface engine
        const meshSurfaceEngine& meshSurface_;

        //- pointer to face patches
        const labelList& facePatch_;

        //- patches at a surface point
        VRWGraph pointPatches_;

        //- labels of corner points
        labelHashSet corners_;

        //- labels of edge points
        labelHashSet edgePoints_;

        //- patch-patches addressing
        List<labelHashSet> patchPatches_;

        //- number of edges attached to a surface point
        labelList nEdgesAtPoint_;

        //- labels of feature edges
        labelHashSet featureEdges_;

    // Private member functions
        //- find corners, edge points and addressing
        void calculateCornersEdgesAndAddressing();

        //- Disallow default bitwise copy construct
        meshSurfacePartitioner(const meshSurfacePartitioner&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfacePartitioner&);

public:

    // Constructors

        //- Construct from meshSurfaceEngine
        meshSurfacePartitioner(const meshSurfaceEngine&);

        //- Construct from meshSurfaceEngine and face patches
        meshSurfacePartitioner
        (
            const meshSurfaceEngine& meshSurface,
            const labelList& facePatch
        );

    // Destructor

        ~meshSurfacePartitioner();

    // Member Functions
        //- return const reference to meshSurfaceEngine
        inline const meshSurfaceEngine& surfaceEngine() const
        {
            return meshSurface_;
        }

        //- return const reference to a list containing boudary patch index
        //- for each boundary face
        inline const labelList& boundaryFacePatches() const
        {
            return facePatch_;
        }

        //- return const reference to a graph of patches attached to a surface
        //- point
        inline const VRWGraph& pointPatches() const
        {
            return pointPatches_;
        }

        //- return labels of corner points (from the list of boundary points)
        inline const labelHashSet& corners() const
        {
            return corners_;
        }

        //- return labels of edge points (from the list of boundary points)
        inline const labelHashSet& edgePoints() const
        {
            return edgePoints_;
        }

        //- return the number of feature edges attached to a boundary point
        inline label numberOfFeatureEdgesAtPoint(const label bpI) const
        {
            return nEdgesAtPoint_[bpI];
        }

        //- return patch-patches addressing
        inline const List<labelHashSet>& patchPatches() const
        {
            return patchPatches_;
        }

        //- return labels of boundary edges which are feature edges
        inline const labelHashSet& featureEdges() const
        {
            return featureEdges_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
