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
    tetFemMatrix<Type>

Description
    Tetrahedral Finite Element matrix.

SourceFiles
    tetFemMatrix.C
    tetFemMatrixSolve.C

\*---------------------------------------------------------------------------*/

#ifndef tetFemMatrix_H
#define tetFemMatrix_H

#include "tetPointFields.H"
#include "elementFields.H"
#include "lduMatrix.H"
#include "tmp.H"
#include "autoPtr.H"
#include "dimensionedTypes.H"
#include "Map.H"
#include "constraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class tetFemMatrix;

template<class Type>
Ostream& operator<<(Ostream&, const tetFemMatrix<Type>&);


/*---------------------------------------------------------------------------*\
                           Class tetFemMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class tetFemMatrix
:
    public refCount,
    public lduMatrix
{
    // Private data

        //- Reference to GeometricField
        GeometricField<Type, tetPolyPatchField, tetPointMesh>& psi_;

        //- Dimension set
        dimensionSet dimensions_;

        //- Source term
        Field<Type> source_;

        //- Are boundary conditions set?
        mutable bool boundaryConditionsSet_;

        //- Equation triangle map
        mutable Map<constraint<Type> > fixedEqns_;

        direction solvingComponent;


    // Private static data

        //- Matrix constraint fill-in
        //  Equals to the estimated fraction of fixed nodes in the matrix
        static const label fixFillIn;

    // Private Member Functions

        //- Add boundary source and diagonal for gradient-type conditions
        void addBoundarySourceDiag();

        //- Store boundary coefficients
        void storeBoundaryCoeffs() const;

        //- Add coupling coefficients
        //  In parallel computations, the consequence of cell-by-cell
        //  matrix assembly is the fact that the "other part" of all the
        //  edge coefficients on the coupled edge and the diagonal
        //  coefficient is assembled on the neighbouring processor.  The
        //  matrix assembly is therefore completed only once the matrix is
        //  completed, by adding the contribution from the neighbouring
        //  domain.
        //  Note: There will be an additional component of
        //  vector-matrix multiply, for the edges originating on the
        //  parallel patch and protruding into the neighbouring domain.
        //  This is treated by doing a "local" multiplication on the
        //  neighbouring side and then adding the result.
        //  HJ, 13/Nov/2001
        void addCouplingCoeffs();

        void eliminateCouplingCoeffs();

        void addCouplingSource(scalarField&) const;

        //- Set component boundary conditions
        void setComponentBoundaryConditions
        (
            const direction,
            scalarField& psiCmpt,
            scalarField& sourceCmpt
        );

        //- Reconstruct matrix
        void reconstructMatrix();

        //- Distribute source.  Given an element field construct a source
        //  by distributing the volume integral over the element points
        tmp<Field<Type> > distributeSource(const Field<Type>&) const;


public:

    // No tetFemSolver: boundary coefficients are not cached
    // HJ, 6/Jul/2007

    ClassName("tetFemMatrix");

    typedef constraint<Type> ConstraintType;


    // Constructors

        //- Construct given a field to solve for
        tetFemMatrix
        (
            GeometricField<Type, tetPolyPatchField, tetPointMesh>&,
            const dimensionSet&
        );

        //- Construct as copy
        tetFemMatrix(const tetFemMatrix<Type>&);

        //- Construct from Istream given field to solve for
        tetFemMatrix
        (
            GeometricField<Type, tetPolyPatchField, tetPointMesh>&,
            Istream&
        );


    // Destructor

        virtual ~tetFemMatrix();


    // Member functions

        // Access

            const GeometricField<Type, tetPolyPatchField, tetPointMesh>&
            psi() const
            {
                return psi_;
            }

            GeometricField<Type, tetPolyPatchField, tetPointMesh>& psi()
            {
                return psi_;
            }

            const dimensionSet& dimensions() const
            {
                return dimensions_;
            }

            Field<Type>& source()
            {
                return source_;
            }

            const Field<Type>& source() const
            {
                return source_;
            }

            //- Does the matrix need a reference level for solution
            //bool needReference();


        // Operations

            //- Set constrained value in a prescribed point
            void addConstraint
            (
                const label vertex,
                const Type& value
            );

            //- Relax matrix (for steady-state solution).
            //  alpha = 1 : diagonally equal
            //  alpha < 1 :    ,,      dominant
            //  alpha = 0 : do nothing
            void relax(const scalar alpha);

            //- Relax matrix (for steady-state solution).
            //  alpha is read from controlDict
            void relax();

            //- Check matrix for diagonal dominance
            void check();

            //- Solve returning the solution statistics.
            //  Convergence tolerance read from dictionary
            lduSolverPerformance solve(const dictionary&);

            //- Solve returning the solution statistics.
            //  Solver controls read from tetFemSolution
            lduSolverPerformance solve();

            //- Return the matrix residual
            tmp<Field<Type> > residual();


    // Member operators

        void operator=(const tetFemMatrix<Type>&);
        void operator=(const tmp<tetFemMatrix<Type> >&);

        void negate();

        void operator+=(const tetFemMatrix<Type>&);
        void operator+=(const tmp<tetFemMatrix<Type> >&);

        void operator-=(const tetFemMatrix<Type>&);
        void operator-=(const tmp<tetFemMatrix<Type> >&);

        void operator+=
        (
            const GeometricField<Type, elementPatchField, elementMesh>&
        );
        void operator+=
        (
            const tmp<GeometricField<Type, elementPatchField, elementMesh> >&
        );

        void operator-=
        (
            const GeometricField<Type, elementPatchField, elementMesh>&
        );

        void operator-=
        (
            const tmp<GeometricField<Type, elementPatchField, elementMesh> >&
        );

        void operator+=(const dimensioned<Type>&);
        void operator-=(const dimensioned<Type>&);

        void operator*=(const dimensioned<scalar>&);

    // Ostream operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const tetFemMatrix<Type>&
        );
};


// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template<class Type>
void checkMethod
(
    const tetFemMatrix<Type>&,
    const tetFemMatrix<Type>&,
    const char*
);


template<class Type>
void checkMethod
(
    const tetFemMatrix<Type>&,
    const GeometricField<Type, elementPatchField, elementMesh>&,
    const char*
);


template<class Type>
void checkMethod
(
    const tetFemMatrix<Type>&,
    const dimensioned<Type>&,
    const char*
);


//- Solve returning the solution statistics given convergence tolerance
//  Solver controls read from Istream
template<class Type>
lduSolverPerformance solve(tetFemMatrix<Type>&, Istream&);


//- Solve returning the solution statistics given convergence tolerance,
//  deleting temporary matrix after solution.
//  Solver controls read from Istream
template<class Type>
lduSolverPerformance solve(const tmp<tetFemMatrix<Type> >&, Istream&);


//- Solve returning the solution statistics given convergence tolerance
//  Convergence tolerance read from controlDict
template<class Type>
lduSolverPerformance solve(tetFemMatrix<Type>&);

//- Solve returning the solution statistics given convergence tolerance,
//  deleting temporary matrix after solution.
//  Convergence tolerance read from controlDict
template<class Type>
lduSolverPerformance solve(const tmp<tetFemMatrix<Type> >&);


// * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * //

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tetFemMatrix<Type>&,
    const tetFemMatrix<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >&,
    const tetFemMatrix<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tetFemMatrix<Type>&,
    const tmp<tetFemMatrix<Type> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >&,
    const tmp<tetFemMatrix<Type> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>&,
    const tetFemMatrix<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >&,
    const tetFemMatrix<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>&,
    const tmp<tetFemMatrix<Type> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >&,
    const tmp<tetFemMatrix<Type> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>&,
    const tetFemMatrix<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >&,
    const tetFemMatrix<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>&,
    const tmp<tetFemMatrix<Type> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >&,
    const tmp<tetFemMatrix<Type> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tetFemMatrix<Type>&,
    const GeometricField<Type, elementPatchField, elementMesh>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >&,
    const GeometricField<Type, elementPatchField, elementMesh>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tetFemMatrix<Type>&,
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >&,
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const GeometricField<Type, elementPatchField, elementMesh>&,
    const tetFemMatrix<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const GeometricField<Type, elementPatchField, elementMesh>&,
    const tmp<tetFemMatrix<Type> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >&,
    const tetFemMatrix<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >&,
    const tmp<tetFemMatrix<Type> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>&,
    const GeometricField<Type, elementPatchField, elementMesh>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >&,
    const GeometricField<Type, elementPatchField, elementMesh>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>&,
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >&,
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const GeometricField<Type, elementPatchField, elementMesh>&,
    const tetFemMatrix<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const GeometricField<Type, elementPatchField, elementMesh>&,
    const tmp<tetFemMatrix<Type> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >&,
    const tetFemMatrix<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >&,
    const tmp<tetFemMatrix<Type> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >&,
    const dimensioned<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const dimensioned<Type>&,
    const tmp<tetFemMatrix<Type> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >&,
    const dimensioned<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const dimensioned<Type>&,
    const tmp<tetFemMatrix<Type> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>&,
    const GeometricField<Type, elementPatchField, elementMesh>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >&,
    const GeometricField<Type, elementPatchField, elementMesh>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>&,
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >&,
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >&,
    const dimensioned<Type>&
);


template<class Type>
tmp<tetFemMatrix<Type> > operator*
(
    const dimensioned<scalar>&,
    const tetFemMatrix<Type>&
);

template<class Type>
tmp<tetFemMatrix<Type> > operator*
(
    const dimensioned<scalar>&,
    const tmp<tetFemMatrix<Type> >&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "tetFemMatrix.C"
#   include "tetFemMatrixSolve.C"
#   include "tetFemMatrixTools.C"
#   include "tetFemMatrixCheck.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
