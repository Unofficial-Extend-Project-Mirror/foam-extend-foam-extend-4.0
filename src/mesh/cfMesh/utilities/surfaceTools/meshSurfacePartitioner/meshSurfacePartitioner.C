/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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
