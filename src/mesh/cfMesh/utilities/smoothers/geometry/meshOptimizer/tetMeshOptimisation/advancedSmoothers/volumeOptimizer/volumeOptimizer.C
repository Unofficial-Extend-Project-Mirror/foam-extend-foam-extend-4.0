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

#include "demandDrivenData.H"
#include "volumeOptimizer.H"
#include "tetrahedron.H"
#include "partTetMeshSimplex.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const vector volumeOptimizer::dirVecs[8] =
    {
        vector(-1.0, -1.0, -1.0),
        vector(1.0, -1.0, -1.0),
        vector(-1.0, 1.0, -1.0),
        vector(1.0, 1.0, -1.0),
        vector(-1.0, -1.0, 1.0),
        vector(1.0, -1.0, 1.0),
        vector(-1.0, 1.0, 1.0),
        vector(1.0, 1.0, 1.0)
    };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volumeOptimizer::volumeOptimizer(partTetMeshSimplex& simplex)
:
    simplexSmoother(simplex)
{}

volumeOptimizer::~volumeOptimizer()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Member functions

void volumeOptimizer::optimizeNodePosition(const scalar tol)
{
    point& p = points_[pointI_];

    if( !bb_.contains(p) )
        p = 0.5 * (bb_.max() + bb_.min());

    const scalar scale = 1.0 / bb_.mag();
    forAll(points_, pI)
        points_[pI] *= scale;
    bb_.min() *= scale;
    bb_.max() *= scale;

    //- find the optimum using divide and conquer
    const scalar func = optimiseDivideAndConquer(tol);
    const point copyP = p;

    //- check if the location can be improved using the steepest descent
    const scalar funcAfter = optimiseSteepestDescent(tol);

    if( funcAfter > func )
        p = copyP;

    //- scale back to the original size
    forAll(points_, pI)
        points_[pI] /= scale;
    bb_.min() /= scale;
    bb_.max() /= scale;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
