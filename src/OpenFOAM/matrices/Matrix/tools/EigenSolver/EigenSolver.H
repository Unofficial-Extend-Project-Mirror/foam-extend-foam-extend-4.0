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
    EigenSolver

Description
    Calculate eigen-values and eigen-vectors of a symmetric dense matrix

SourceFiles
    EigenSolver.C

\*---------------------------------------------------------------------------*/

#ifndef EigenSolver_H
#define EigenSolver_H

#include "SquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class EigenSolver Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class EigenSolver
{
    // Private data

        //- Eigenvalues
        List<T> values_;

        //- Eigenvectors
        List<List<T> > vectors_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        EigenSolver(const EigenSolver&);

        //- Disallow default bitwise assignment
        void operator=(const EigenSolver&);


        //- Check matrix for shape and symmetry
        void checkMatrix(const SquareMatrix<T>& mtx) const;

        //- Factorise into eigen-values and eigen-vectors
        void factorise(const SquareMatrix<T>& mtx);

        //- Rotate the matrix
        inline void rotate
        (
            SquareMatrix<T>& a,
            const T s,
            const T tau,
            const label i,
            const label j,
            const label k,
            const label l
        ) const;


public:

    // Static data members


    // Constructors

        //- Construct from matrix
        EigenSolver(const SquareMatrix<T>& mtx);


    // Destructor - default


    // Member Functions

        //- Return nth eigen value
        T eigenValue(const label n) const;

        // Return nth eigen vector
        const List<T>& eigenVector(const label n) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EigenSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
