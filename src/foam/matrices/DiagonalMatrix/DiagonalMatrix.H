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
    Foam::DiagonalMatrix<Type>

Description
    DiagonalMatrix<Type> is a 2D diagonal matrix of objects
    of type Type, size nxn

SourceFiles
    DiagonalMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef DiagonalMatrix_H
#define DiagonalMatrix_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * *  * * * * * * Class Forward declaration  * * * * * * * * * * * //

template<class Form, class Type> class Matrix;

/*---------------------------------------------------------------------------*\
                           Class DiagonalMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DiagonalMatrix
:
    public List<Type>
{
public:

    // Constructors

        //- Construct from diagonal component of a Matrix
        template<class Form>
        DiagonalMatrix<Type>(const Matrix<Form, Type>&);

        //- Construct empty from size
        DiagonalMatrix<Type>(const label size);

        //- Construct from size and a value
        DiagonalMatrix<Type>(const label, const Type&);


    // Member functions

        //- Invert the diaganol matrix and return itself
        DiagonalMatrix<Type>& invert();
};


// Global functions

//- Return the diagonal Matrix inverse
template<class Type>
DiagonalMatrix<Type> inv(const DiagonalMatrix<Type>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DiagonalMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
