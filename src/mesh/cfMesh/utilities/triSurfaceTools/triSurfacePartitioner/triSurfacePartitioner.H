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
    triSurfacePartitioner

Description
    Provides information regarding surface partitions on the surface
    triangulation, and connectivity between various surface partitions.

SourceFiles
    triSurfacePartitioner.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfacePartitioner_H
#define triSurfacePartitioner_H

#include "triSurf.H"
#include "DynList.H"
#include "HashSet.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class triSurfacePartitioner Declaration
\*---------------------------------------------------------------------------*/

class triSurfacePartitioner
{
    // Private data
        //- reference to triSurf
        const triSurf& surface_;

        //- corner nodes
        labelList corners_;
        List<DynList<label> > cornerPatches_;

        //- information which partitions share an edge with a given partition
        List<labelHashSet> patchPatches_;

        //- edge partitions
        labelList edgeGroups_;

        //- information which edge groups share a corner
        //- with a given group
        List<labelHashSet> edgeGroupEdgeGroups_;

        //- edge groups between surface patches
        std::map<std::pair<label, label>, labelHashSet> patchesEdgeGroups_;

        //- corners shared by edge groups
        std::map<std::pair<label, label>, labelHashSet> edgeGroupsCorners_;

    // Private member functions
        //- calculate patch addressing
        void calculatePatchAddressing();

        //- find surface corners
        void calculateCornersAndAddressing();

        //- calculate patch-patches addressing
        void calculatePatchPatches();

        //- calculate edge groups
        void calculateEdgeGroups();

        //- calculate surface patch to edge groups addressing
        void calculatePatchToEdgeGroups();

        //- calculate edge groups to corner addressing
        void calculateEdgeGroupsToCorners();

        //- Disallow default bitwise copy construct
        triSurfacePartitioner(const triSurfacePartitioner&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfacePartitioner&);

public:

    // Constructors

        //- Construct from triSurf
        triSurfacePartitioner(const triSurf& surface);

    // Destructor

        ~triSurfacePartitioner();


    // Member Functions
        //- return corner nodes
        const labelList& corners() const;

        //- return corner patches
        const List<DynList<label> >& cornerPatches() const;

        //- return patch-patches addressing
        const List<labelHashSet>& patchPatches() const;

        //- return edge groups. Edges which are not
        //- feature edges are set to -1
        const labelList& edgeGroups() const;

        //- Edge group - edge groups addressing
        const List<labelHashSet>& edgeGroupEdgeGroups() const;

        //- Return groups of feature edges shared by the given patches
        void edgeGroupsSharedByPatches
        (
            const label partition1,
            const label partition2,
            DynList<label>& edgePartitions
        ) const;

        //- return corners shared shared by the given edge groups
        void cornersSharedByEdgeGroups
        (
            const label edgePartition1,
            const label edgePartition2,
            DynList<label>& corners
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
