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
    matrix2D

Description
    Implementation of 2 x 2 matrix

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef matrix2D_H
#define matrix2D_H

#include "scalar.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class matrix2D Declaration
\*---------------------------------------------------------------------------*/

class matrix2D
: public FixedList<FixedList<scalar, 2>, 2>
{
    // Private members
        scalar det_;
        bool calculatedDet_;

    // Private member functions
        //- calculate matrix determinant
        inline void calculateDeterminant();

public:

    // Constructors
        //- Null constructor
        explicit inline matrix2D();

        //- Copy constructor
        inline matrix2D(const matrix2D&);

    // Destructor
        inline ~matrix2D();

    // Member functions
        //- return matrix determinant
        inline scalar determinant();

        //- return inverse matrix
        inline matrix2D inverse();

        //- find the solution of the system with the given rhs
        inline FixedList<scalar, 2> solve
        (
            const FixedList<scalar, 2>& source
        );

        //- return the first component of the solution vector
        inline scalar solveFirst(const FixedList<scalar, 2>& source);

        //- return the second component of the solution vector
        inline scalar solveSecond(const FixedList<scalar, 2>& source);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "matrix2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
