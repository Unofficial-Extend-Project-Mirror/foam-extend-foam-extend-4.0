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
    Foam::simpleMatrix

Description
    A simple square matrix solver with scalar coefficients.
    Simple dense matrix class with direct solvers.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    simpleMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef simpleMatrix_H
#define simpleMatrix_H

#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class simpleMatrix;

template<class Type>
simpleMatrix<Type> operator+
(
    const simpleMatrix<Type>&,
    const simpleMatrix<Type>&
);

template<class Type>
simpleMatrix<Type> operator-
(
    const simpleMatrix<Type>&,
    const simpleMatrix<Type>&
);

template<class Type>
simpleMatrix<Type> operator*
(
    const scalar,
    const simpleMatrix<Type>&
);

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const simpleMatrix<Type>&
);


/*---------------------------------------------------------------------------*\
                           Class simpleMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class simpleMatrix
:
    public scalarSquareMatrix
{
    // Private data

        Field<Type> source_;


public:

    // Constructors

        //- Construct given size
        //  Note: this does not initialise the coefficients or the source.
        simpleMatrix(const label);

        //- Construct given size and initial values for the
        //  coefficients and source
        simpleMatrix(const label, const scalar, const Type&);

        //- Construct from components
        simpleMatrix(const scalarSquareMatrix&, const Field<Type>&);

        //- Construct from Istream
        simpleMatrix(Istream&);

        //- Construct as copy
        simpleMatrix(const simpleMatrix<Type>&);


    // Member Functions

        // Access

            //- Return access to the source
            Field<Type>& source()
            {
                return source_;
            }

            //- Return const-access to the source
            const Field<Type>& source() const
            {
                return source_;
            }


        //- Solve the matrix using Gaussian elimination with pivoting
        //  and return the solution
        Field<Type> solve() const;

        //- Solve the matrix using LU decomposition with pivoting
        //  and return the solution
        Field<Type> LUsolve() const;


    // Member Operators

        void operator=(const simpleMatrix<Type>&);


    // Ostream Operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const simpleMatrix<Type>&
        );
};


// Global operators

template<class Type>
simpleMatrix<Type> operator+
(
    const simpleMatrix<Type>&,
    const simpleMatrix<Type>&
);

template<class Type>
simpleMatrix<Type> operator-
(
    const simpleMatrix<Type>&,
    const simpleMatrix<Type>&
);

template<class Type>
simpleMatrix<Type> operator*
(
    const scalar,
    const simpleMatrix<Type>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "simpleMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
