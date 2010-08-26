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
    Foam::Matrix

Description
    A templated 2D matrix of objects of \<T\>, where the n x m matrix
    dimensions are known and used for subscript bounds checking, etc.

SourceFiles
    Matrix.C
    MatrixI.H
    MatrixIO.C

\*---------------------------------------------------------------------------*/

#ifndef Matrix_H
#define Matrix_H

#include "List.H"
#include "label.H"
#include "bool.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T> class Matrix;

template<class T> const T& max(const Matrix<T>&);
template<class T> const T& min(const Matrix<T>&);

template<class T> Matrix<T> operator-(const Matrix<T>&);
template<class T> Matrix<T> operator+(const Matrix<T>&, const Matrix<T>&);
template<class T> Matrix<T> operator-(const Matrix<T>&, const Matrix<T>&);
template<class T> Matrix<T> operator*(const scalar, const Matrix<T>&);

template<class T> Istream& operator>>(Istream&, Matrix<T>&);
template<class T> Ostream& operator<<(Ostream&, const Matrix<T>&);


/*---------------------------------------------------------------------------*\
                           Class Matrix Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Matrix
{
    // Private data

        //- Row pointers
        T** __restrict__ v_;

        //- Number of rows and columns in Matrix.
        label n_, m_;

        //- Allocate the storage for the row-pointers and the data
        //  and set the row pointers
        void allocate();


public:

    // Constructors

        //- Null constructor.
        inline Matrix();

        //- Construct given number of rows and columns.
        Matrix(const label n, const label m);

        //- Construct with given number of rows and columns
        //  and value for all elements.
        Matrix(const label n, const label m, const T&);

        //- Copy constructor.
        Matrix(const Matrix<T>&);

        //- Construct from Istream.
        Matrix(Istream&);

        //- Clone
        inline autoPtr<Matrix<T> > clone() const;


    // Destructor

        ~Matrix();


    // Member functions

        //- Return a null Matrix
        static const Matrix<T>& null();


        // Access

            //- Return the number of rows
            inline label n() const;

            //- Return the number of columns
            inline label m() const;

            //- Return the number of elements in matrix (n*m)
            inline label size() const;


        // Check

            //- Check index i is within valid range (0 ... n-1).
            inline void checki(const label i) const;

            //- Check index j is within valid range (0 ... m-1).
            inline void checkj(const label j) const;


        // Edit

            //- Clear the Matrix, i.e. set sizes to zero.
            void clear();

            //- Transfer the contents of the argument Matrix into this Matrix
            //  and annull the argument Matrix.
            void transfer(Matrix<T>&);


    // Member operators

        //- Return subscript-checked element of Matrix.
        inline T* operator[](const label);

        //- Return subscript-checked element of constant Matrix.
        inline const T* operator[](const label) const;

        //- Assignment operator. Takes linear time.
        void operator=(const Matrix<T>&);

        //- Assignment of all entries to the given value
        void operator=(const T&);


    // Friend functions

        friend const T& max<T>(const Matrix<T>&);
        friend const T& min<T>(const Matrix<T>&);


    // Friend operators

        friend Matrix<T> operator-<T>(const Matrix<T>&);
        friend Matrix<T> operator+<T>(const Matrix<T>&, const Matrix<T>&);
        friend Matrix<T> operator-<T>(const Matrix<T>&, const Matrix<T>&);
        friend Matrix<T> operator*<T>(const scalar, const Matrix<T>&);


    // IOstream operators

        //- Read Matrix from Istream, discarding contents of existing Matrix.
        friend Istream& operator>> <T>(Istream&, Matrix<T>&);

        // Write Matrix to Ostream.
        friend Ostream& operator<< <T>(Ostream&, const Matrix<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "MatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Matrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
