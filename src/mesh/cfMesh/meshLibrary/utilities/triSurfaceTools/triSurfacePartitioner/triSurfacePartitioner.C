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

Description

\*---------------------------------------------------------------------------*/

#include "triSurfacePartitioner.H"
#include "demandDrivenData.H"

# ifdef DEBUGPartitioner
#include <sstream>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfacePartitioner::triSurfacePartitioner
(
    const triSurf& surface
)
:
    surface_(surface),
    corners_(),
    cornerPatches_(),
    patchPatches_(surface.patches().size()),
    edgeGroups_(),
    edgeGroupEdgeGroups_(),
    patchesEdgeGroups_(),
    edgeGroupsCorners_()
{
    calculatePatchAddressing();
}

triSurfacePartitioner::~triSurfacePartitioner()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const labelList& triSurfacePartitioner::corners() const
{
    return corners_;
}

const List<DynList<label> >& triSurfacePartitioner::cornerPatches() const
{
    return cornerPatches_;
}

const List<labelHashSet>& triSurfacePartitioner::patchPatches() const
{
    return patchPatches_;
}

const labelList& triSurfacePartitioner::edgeGroups() const
{
    return edgeGroups_;
}

const List<labelHashSet>& triSurfacePartitioner::edgeGroupEdgeGroups() const
{
    return edgeGroupEdgeGroups_;
}

void triSurfacePartitioner::edgeGroupsSharedByPatches
(
    const label patch1,
    const label patch2,
    DynList<label>& edgeGroups
) const
{
    edgeGroups.clear();

    std::pair<label, label> pp
    (
        Foam::min(patch1, patch2),
        Foam::max(patch1, patch2)
    );

    std::map<std::pair<label, label>, labelHashSet>::const_iterator it =
        patchesEdgeGroups_.find(pp);

    if( it != patchesEdgeGroups_.end() )
    {
        const labelHashSet& eGroups = it->second;

        forAllConstIter(labelHashSet, eGroups, it)
            edgeGroups.append(it.key());
    }
}

void triSurfacePartitioner::cornersSharedByEdgeGroups
(
    const label edgeGroup1,
    const label edgeGroup2,
    DynList<label>& corners
) const
{
    corners.clear();

    std::pair<label, label> ep
    (
        Foam::min(edgeGroup1, edgeGroup2),
        Foam::max(edgeGroup1, edgeGroup2)
    );

    std::map<std::pair<label, label>, labelHashSet>::const_iterator it =
        edgeGroupsCorners_.find(ep);

    if( it != edgeGroupsCorners_.end() )
    {
        const labelHashSet& corn = it->second;

        forAllConstIter(labelHashSet, corn, it)
            corners.append(it.key());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
