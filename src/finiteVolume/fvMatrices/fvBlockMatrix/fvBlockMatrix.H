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
    Foam::fvBlockMatrix

Description
    fvBlockMatrix is an extension of fvMatrix for block coupled types. It holds
    general insertion and retrieval tools for block coupling and specific
    functions for p-U coupling.

Author
    Vuko Vukcevic, FMENA Zagreb.
    Update by Hrvoje Jasak

SourceFiles
    fvBlockMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef fvBlockMatrix_H
#define fvBlockMatrix_H

#include "BlockLduSystem.H"
#include "fvMatrices.H"
#include "blockLduSolvers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class fvBlockMatrix;

template<class Type>
Ostream& operator<<(Ostream&, const fvBlockMatrix<Type>&);


/*---------------------------------------------------------------------------*\
                        Class fvBlockMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvBlockMatrix
:
    public BlockLduSystem<Type, Type>
{
    // Private data

        //- Internal field reference
        GeometricField<Type, fvPatchField, volMesh>& psi_;


    // Private member functions

        // Insertion functions for fvMatrix<matrixType> into diagonal positions

            //- Insert internal field into this fvBlockMatrix
            template<class fieldType>
            void insertSolutionVector
            (
                const direction dir,
                const Field<fieldType>& xSingle
            );

            //- Insert matrix diagonal and source into this fvBlockMatrix
            template<class matrixType>
            void insertDiagSource
            (
                const direction dir,
                fvMatrix<matrixType>& matrix
            );

            //- Insert upper and lower part into this fvBlockMatrix
            template<class matrixType>
            void insertUpperLower
            (
                const direction dir,
                const fvMatrix<matrixType>& matrix
            );

            //- Update coupling coefficients in this fvBlockMatrix
            template<class matrixType>
            void updateCouplingCoeffs
            (
                const direction dir,
                const fvMatrix<matrixType>& matrix
            );


        // Insertion functions for fvScalarMatrix into off-diagonal positions
        // (coupling matrices)

            //- Insert coupling matrix diag element into this fvBlockMatrix
            void insertCouplingDiag
            (
                const direction dirI,
                const direction dirJ,
                const scalarField& coeffIJ
            );

            //- Insert coupling matrix diag and source into this fvBlockMatrix
            void insertCouplingDiagSource
            (
                const direction dirI,
                const direction dirJ,
                const fvScalarMatrix& matrix
            );

            //- Insert coupling matrix lower and upper into this fvBlockMatrix
            void insertCouplingUpperLower
            (
                const direction dirI,
                const direction dirJ,
                const fvScalarMatrix& matrix
            );


        // Pressure - velocity (p-U) coupling specific functions

            //- Insert BlockLduSystem (obtained by implicit grad/div operator)
            // into this fvBlockMatrix
            template<class blockType, class fieldType>
            void insertBlock
            (
                const direction dirI,
                const direction dirJ,
                const BlockLduSystem<blockType, fieldType>& blockSystem,
                const bool incrementColumnDir
            );

            //- Insert source and coupling coeffs of BlockLduSystem
            // (eg. obtained by implicit grad/div operator)
            template<class blockType, class fieldType>
            void insertBoundaryContributions
            (
                const direction dirI,
                const direction dirJ,
                const BlockLduSystem<blockType, fieldType>& blockSystem,
                const bool incrementColumnDir
            );


public:

    // Constructors

        //- Construct given a field to solve for
        fvBlockMatrix(GeometricField<Type, fvPatchField, volMesh>&);

        //- Construct as copy
        fvBlockMatrix(const fvBlockMatrix<Type>&);


    // Destructor

       virtual ~fvBlockMatrix()
       {}


    // Member functions

        // Access

            //- Access to GeometricField
            GeometricField<Type, fvPatchField, volMesh>& psi()
            {
                return psi_;
            }

            //- const reference to GeometricField
            const GeometricField<Type, fvPatchField, volMesh>& psi() const
            {
                return psi_;
            }


        // Insertion and retrieval public tools

            //- Retrieve part of internal field from this fvBlockMatrix
            template<class fieldType>
            void retrieveSolution
            (
                const direction dir,
                Field<fieldType>& xSingle
            ) const;

            //- Insert matrix into this fvBlockMatrix
            template<class matrixType>
            void insertEquation
            (
                const direction dir,
                fvMatrix<matrixType>& matrix
            );

            //- Insert existing block system (obtained by implicit grad/div
            // operator) into this fvBlockMatrix
            template<class blockType, class fieldType>
            void insertBlockCoupling
            (
                const direction dirI,
                const direction dirJ,
                const BlockLduSystem<blockType, fieldType>& blockSystem,
                const bool incrementColumnDir
            );

            //- Insert diagonal only equation coupling into this fvBlockMatrix
            //  Source compensation is done in function updateSourceCoupling()
            //  after all coupling terms are added.  HJ, 27/Apr/2015
            void insertEquationCoupling
            (
                const direction dirI,
                const direction dirJ,
                const scalarField& coeffIJ
            );

            //- Insert scalar equation coupling into this fvBlockMatrix
            void insertEquationCoupling
            (
                const direction dirI,
                const direction dirJ,
                const fvScalarMatrix& matrix
            );

            //- Add field into block field
            void blockAdd
            (
                const direction dir,
                const scalarField& xSingle,
                Field<Type>& blockX
            );

            //- Update coupling of block system.
            //  Subtracts the block-coefficient coupling as specified by the
            //  user from the source, leaving the implicit update given by
            //  linearisation
            void updateSourceCoupling();

            //- Insert adjoint convection tensor term. Hard coded to correspond
            //  to Gauss gradient discretisation with linear interpolation. This
            //  needs to be reorganised. VV, 7/Apr/2016.
            void insertAdjointConvection
            (
                const direction UEqnDir,
                const volVectorField& U,
                const volVectorField& UStar
            );

            //- Insert Picard tensor term that comes from Picard linearisation
            // of convection term in momentum equation. VV, 21/July/2014.
            void insertPicardTensor
            (
                const direction UEqnDir,
                const volVectorField& U,
                const surfaceScalarField& phi
            );


        // Solver calls for fvBlockMatrix

            //- Solve returning the solution statistics.
            // Use the given solver controls
            BlockSolverPerformance<Type> solve(const dictionary&);

            //- Solve returning the solution statistics.
            // Solver controls read from fvSolution
            BlockSolverPerformance<Type> solve();


    // Member operators

        void operator=(const fvBlockMatrix<Type>&);

        void negate();

        void operator+=(const fvBlockMatrix<Type>&);
        void operator-=(const fvBlockMatrix<Type>&);

        void operator*=(const scalarField&);
        void operator*=(const scalar);


    // Ostream operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const fvBlockMatrix<Type>&
        );
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvBlockMatrix.C"
#   include "fvBlockMatrixOperations.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
