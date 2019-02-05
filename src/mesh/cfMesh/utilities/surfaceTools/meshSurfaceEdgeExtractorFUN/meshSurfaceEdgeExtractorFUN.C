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

#include "meshSurfaceEdgeExtractorFUN.H"
#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

meshSurfaceEngine& meshSurfaceEdgeExtractorFUN::surfaceEngine()
{
    # ifdef USE_OMP
    if( omp_in_parallel() )
        FatalErrorIn
        (
            "meshSurfaceEngine& meshSurfaceEdgeExtractorFUN::surfaceEngine()"
        ) << "Cannot create surface engine with a parallel region"
            << exit(FatalError);
    # endif

    if( !surfaceEnginePtr_ )
        surfaceEnginePtr_ = new meshSurfaceEngine(mesh_);

    return *surfaceEnginePtr_;
}

void meshSurfaceEdgeExtractorFUN::clearOut()
{
    deleteDemandDrivenData(surfaceEnginePtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh and octree
meshSurfaceEdgeExtractorFUN::meshSurfaceEdgeExtractorFUN
(
    polyMeshGen& mesh,
    const meshOctree& octree,
    const bool createWrapperSheet
)
:
    mesh_(mesh),
    meshOctree_(octree),
    surfaceEnginePtr_(nullptr),
    createWrapperSheet_(createWrapperSheet)
{
    if( Pstream::parRun() )
        FatalErrorIn
        (
            "meshSurfaceEdgeExtractorFUN::meshSurfaceEdgeExtractorFUN"
            "(polyMeshGen&, const meshOctree&)"
        ) << "Cannot run in parallel!" << exit(FatalError);

    createBasicFundamentalSheets();

    smoothMeshSurface();

    remapBoundaryPoints();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceEdgeExtractorFUN::~meshSurfaceEdgeExtractorFUN()
{
    deleteDemandDrivenData(surfaceEnginePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
