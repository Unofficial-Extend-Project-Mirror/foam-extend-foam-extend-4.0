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

#include "meshSurfacePartitioner.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfacePartitioner::meshSurfacePartitioner
(
    const meshSurfaceEngine& meshSurface
)
:
    meshSurface_(meshSurface),
    facePatch_(meshSurface.boundaryFacePatches()),
    pointPatches_(),
    corners_(),
    edgePoints_(),
    patchPatches_(),
    nEdgesAtPoint_(),
    featureEdges_()
{
    calculateCornersEdgesAndAddressing();
}

meshSurfacePartitioner::meshSurfacePartitioner
(
    const meshSurfaceEngine& meshSurface,
    const labelList& facePatch
)
:
    meshSurface_(meshSurface),
    facePatch_(facePatch),
    pointPatches_(),
    corners_(),
    edgePoints_(),
    patchPatches_(),
    nEdgesAtPoint_(),
    featureEdges_()
{
    calculateCornersEdgesAndAddressing();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfacePartitioner::~meshSurfacePartitioner()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
