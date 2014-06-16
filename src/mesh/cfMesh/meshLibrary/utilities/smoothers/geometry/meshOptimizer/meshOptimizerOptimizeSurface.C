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

#include "demandDrivenData.H"
#include "meshOptimizer.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceMapper.H"
#include "meshOctree.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOptimizer::optimizeSurface(const meshOctree& octree)
{
    Info << "Optimizing positions of surface nodes" << endl;
    
    meshSurfaceEngine& mse = const_cast<meshSurfaceEngine&>(meshSurface());
    meshSurfaceOptimizer surfaceOptimizer(mse, octree);

    surfaceOptimizer.optimizeSurface();
    
    meshSurfaceMapper(mse, octree).mapVerticesOntoSurfacePatches();
    
    clearSurface();

    Info << "Finished optimizing positions of surface nodes" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
