/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::RectangularMatrix

Description
    A templated 2D rectangular matrix of objects of \<T\>, where the
    n x m matrix dimension is known and used for subscript bounds checking etc.

SourceFiles
    RectangularMatrixI.H
    RectangularMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef RectangularMatrix_H
#define RectangularMatrix_H

#include "Matrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Matrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class RectangularMatrix
:
    public Matrix<RectangularMatrix<Type>, Type>
{

public:

    // Constructors

        //- Null constructor.
        inline RectangularMatrix();

        //- Construct given number of rows and columns,
        inline RectangularMatrix(const label m, const label n);

        //- Construct with given number of rows and columns
        //  and value for all elements.
        inline RectangularMatrix(const label m, const label n, const Type&);

        //- Construct from Istream.
        inline RectangularMatrix(Istream&);

        //- Clone
        inline autoPtr<RectangularMatrix<Type> > clone() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "RectangularMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
