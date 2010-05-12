/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::simpleMatrix

Description
    Simple dense matrix class with direct solvers.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    simpleMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef simpleMatrix_H
#define simpleMatrix_H

#include "scalarMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T>
class simpleMatrix;

template<class T>
simpleMatrix<T> operator+
(
    const simpleMatrix<T>&,
    const simpleMatrix<T>&
);

template<class T>
simpleMatrix<T> operator-
(
    const simpleMatrix<T>&,
    const simpleMatrix<T>&
);

template<class T>
simpleMatrix<T> operator*
(
    const scalar,
    const simpleMatrix<T>&
);

template<class T>
Ostream& operator<<
(
    Ostream&,
    const simpleMatrix<T>&
);


/*---------------------------------------------------------------------------*\
                           Class simpleMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class simpleMatrix
:
    public scalarMatrix
{
    // Private data

        Field<T> source_;


public:

    // Constructors

        //- Construct given size
        simpleMatrix(const label);

        //- Construct from components
        simpleMatrix(const scalarMatrix&, const Field<T>&);

        //- Construct from Istream
        simpleMatrix(Istream&);

        //- Construct as copy
        simpleMatrix(const simpleMatrix<T>&);


    // Member Functions

        // Access

            Field<T>& source()
            {
                return source_;
            }

            const Field<T>& source() const
            {
                return source_;
            }


        //- Solve the matrix using Gaussian elimination with pivoting
        //  and return the solution
        Field<T> solve() const;

        //- Solve the matrix using LU decomposition with pivoting
        //  and return the solution
        Field<T> LUsolve() const;


    // Member Operators

        void operator=(const simpleMatrix<T>&);


    // Friend Operators

        friend simpleMatrix<T> operator+ <T>
        (
            const simpleMatrix<T>&,
            const simpleMatrix<T>&
        );

        friend simpleMatrix<T> operator- <T>
        (
            const simpleMatrix<T>&,
            const simpleMatrix<T>&
        );

        friend simpleMatrix<T> operator* <T>
        (
            const scalar,
            const simpleMatrix<T>&
        );


    // Ostream Operator

        friend Ostream& operator<< <T>
        (
            Ostream&,
            const simpleMatrix<T>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "simpleMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
