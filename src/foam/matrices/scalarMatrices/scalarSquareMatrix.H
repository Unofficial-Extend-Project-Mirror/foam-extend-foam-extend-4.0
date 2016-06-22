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
    Foam::scalarSquareMatrix

Description
    Foam::scalarSquareMatrix

SourceFiles
    scalarSquareMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef scalarSquareMatrix_H
#define scalarSquareMatrix_H

#include "SquareMatrix.H"
#include "scalarField.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class scalarSquareMatrix Declaration
\*---------------------------------------------------------------------------*/

class scalarSquareMatrix
:
    public SquareMatrix<scalar>
{

public:

    // Constructors

        //- Construct null
        scalarSquareMatrix();

        //- Construct given size
        scalarSquareMatrix(const label);

        //- Construct given size and coefficient value
        scalarSquareMatrix(const label, const scalar v);

        //- Construct as copy
        scalarSquareMatrix(const scalarSquareMatrix&);

        //- Construct from Istream
        scalarSquareMatrix(Istream&);


    // Member Functions

        //- Solve the matrix using Gaussian elimination with pivoting,
        //  returning the solution in the source.
        //  Matrix given as a reference is destroyed during solution
        template<class T>
        static void solve(scalarSquareMatrix& matrix, Field<T>& source);

        //- Solve the matrix using Gaussian elimination with pivoting
        //  and return the solution
        template<class Type>
        void solve(Field<Type>& psi, const Field<Type>& source) const;


        //- LU decompose the matrix with pivoting
        static void LUDecompose
        (
            scalarSquareMatrix& matrix,
            labelList& pivotIndices
        );

        //- LU back-substitution with given source, returning the solution
        //  in the source
        template<class Type>
        static void LUBacksubstitute
        (
            const scalarSquareMatrix& luMmatrix,
            const labelList& pivotIndices,
            Field<Type>& source
        );

        //- Invert the matrix using LU decomposition with pivoting
        scalarSquareMatrix LUinvert() const;

        //- Solve the matrix using LU decomposition with pivoting
        //  returning the LU form of the matrix and the solution in the source
        template<class Type>
        static void LUsolve(scalarSquareMatrix& matrix, Field<Type>& source);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "scalarSquareMatrixTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
