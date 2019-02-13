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

Class
    volumeOptimizer

Description
    Mesh smoothing without any topological changes. The vertex is put into
    location which minimises the objective function.

SourceFiles
    volumeOptimizer.C

\*---------------------------------------------------------------------------*/

#ifndef volumeOptimizer_H
#define volumeOptimizer_H

#include "simplexSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class partTetMeshSimplex;

/*---------------------------------------------------------------------------*\
                    Class volumeOptimizer Declaration
\*---------------------------------------------------------------------------*/

//- class for volume optimizer
class volumeOptimizer
:
    public simplexSmoother
{
    // Private static data
        //- direction vectors for divide and conquer algorithm
        static const vector dirVecs[8];

    // Private member functions
        //- evaluate functional
        scalar evaluateFunc() const;

        //- find appropriate value of K
        scalar evaluateStabilisationFactor() const;

        //- evaluate gradients of the functional
        void evaluateGradientsExact(vector&, tensor&) const;

        //- optimize position using a divide and conquer algorithm
        scalar optimiseDivideAndConquer(const scalar tol);

        //- optimise using the steepest descent
        scalar optimiseSteepestDescent(const scalar tol);

    public:

    // Constructor
        //- construct from partTetMeshSimplex
        volumeOptimizer(partTetMeshSimplex& simplex);

    // Destructor
        ~volumeOptimizer();

    // Member functions
        //- find the best position for the node
        void optimizeNodePosition(const scalar tol = 0.001);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
