/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    BlockLduMatrix

Description
    BlockLduMatrix is a general matrix class in which the coefficients are
    stored as three arrays, one for the upper triangle, one for the
    lower triangle and a third for the diagonal.  Addressing object must
    be supplied for the upper and lower triangles.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    BlockLduMatrix.C
    BlockLduMatrixOperations.C
    BlockLduMatrixUpdateInterfaces.C
    BlockLduMatrixATmul.C
    BlockLduMatrixHOps.C
    BlockLduMatrixDecouple.C
    BlockLduMatrixDecoupledHOps.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduMatrix_H
#define BlockLduMatrix_H

#include "coeffFields.H"
#include "FieldField.H"
#include "lduMesh.H"
// #include "BlockLduInterfaceField.H"
#include "BlockLduInterfaceFieldPtrsList.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class BlockLduMatrix;

template<class Type>
Ostream& operator<<(Ostream&, const BlockLduMatrix<Type>&);

template<class Type>
class BlockConstraint;


/*---------------------------------------------------------------------------*\
                        Class BlockLduMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockLduMatrix
:
    public refCount
{
public:

    // Public data types

        typedef CoeffField<Type> TypeCoeffField;
        typedef Field<Type> TypeField;
        typedef BlockConstraint<Type> ConstraintType;


private:

    // Private data

        // LDU mesh reference
        const lduMesh& lduMesh_;


        // Block matrix elements

            //- Diagonal coefficients
            CoeffField<Type>* diagPtr_;

            //- Upper triangle coefficients.  Also used for symmetric matrix
            CoeffField<Type>* upperPtr_;

            //- Lower triangle coefficients
            CoeffField<Type> *lowerPtr_;


        // Coupling

            //- List of coupled interfaces
            typename BlockLduInterfaceFieldPtrsList<Type>::Type interfaces_;

            //- Coupled interface coefficients, upper
            FieldField<CoeffField, Type> coupleUpper_;

            //- Coupled interface coefficients, lower
            FieldField<CoeffField, Type> coupleLower_;


        // Constraints

            //- Equation triangle map
            mutable Map<ConstraintType> fixedEqns_;


    // Private static data

        //- Matrix constraint fill-in
        //  Equals to the estimated fraction of fixed nodes in the matrix
        static const label fixFillIn;


    // Private member functions

        // Decoupled versions of nmatrix operations

            //- Sum off-diagonal coefficients and add to diagonal,
            //  decoupled version
            void decoupledSumDiag();

            //- Sum negative off-diagonal coefficients and add to diagonal,
            //  decoupled version
            void decoupledNegSumDiag();

            //- Check matrix for diagonal dominance, decoupled version
            void decoupledCheck() const;

            //- Relax matrix, decoupled version
            void decoupledRelax
            (
                const TypeField& x,
                TypeField& b,
                const scalar alpha
            );

            //- Matrix scaling with scalar field, decoupled version
            void decoupledMultEqOp(const scalarField& sf);

            //- Matrix multiplication without coupled interface update,
            //  Decoupled version
            void decoupledAmulCore
            (
                TypeField& Ax,
                const TypeField& x
            ) const;

            //- Matrix transpose multiplication without coupled
            //  interface update
            //  Decoupled version
            void decoupledTmulCore
            (
                TypeField& Tx,
                const TypeField& x
            ) const;

            //- Return L-U vector-matrix multiplication in row-form,
            //  Decoupled version
            tmp<TypeField> decoupledH(const TypeField& x) const;

            //- Return L-U  vector-matrix multiplication in off-diagonal form,
            //  Decoupled version
            tmp<TypeField> decoupledFaceH(const TypeField& x) const;


protected:

    // Access to constraints

        //- Return constraint map
        const Map<ConstraintType>& fixedEqns() const
        {
            return fixedEqns_;
        }


        //- Return access constraint map
        Map<ConstraintType>& fixedEqns()
        {
            return fixedEqns_;
        }


public:

    //- Runtime type information
    TypeName("BlockLduMatrix");


    // Constructors

        //- Construct given addressing
        explicit BlockLduMatrix(const lduMesh&);

        //- Construct as copy
        BlockLduMatrix(const BlockLduMatrix<Type>&);

        //- Construct as copy or re-use as specified.
        BlockLduMatrix(BlockLduMatrix<Type>&, bool reUse);


    // Destructor

        virtual ~BlockLduMatrix();


    // Member functions

        // Access to addressing

            //- Return the LDU mesh from which the addressing is obtained
            const lduMesh& mesh() const
            {
                return lduMesh_;
            }

            //- Return the LDU addressing
            const lduAddressing& lduAddr() const
            {
                return lduMesh_.lduAddr();
            }

            //- Return the patch evaluation schedule
            const lduSchedule& patchSchedule() const
            {
                return lduAddr().patchSchedule();
            }


        // Access to coefficients

            //- Return access to diagonal coefficients
            TypeCoeffField& diag();

            //- Return diagonal coefficients
            const TypeCoeffField& diag() const;

            //- Return access to upper coefficients
            //  Also used for symmetric matrices
            TypeCoeffField& upper();

            //- Return upper coefficients
            //  Also used for symmetric matrices
            const TypeCoeffField& upper() const;

            //- Return access to lower coefficients
            TypeCoeffField& lower();

            //- Return lower coefficients
            const TypeCoeffField& lower() const;

            //- Return access to coupled interface coefficients, upper
            FieldField<CoeffField, Type>& coupleUpper()
            {
                return coupleUpper_;
            }

            //- Return coupled interface coefficients, upper
            const FieldField<CoeffField, Type>& coupleUpper() const
            {
                return coupleUpper_;
            }

            //- Return access to coupled interface coefficients, lower
            FieldField<CoeffField, Type>& coupleLower()
            {
                return coupleLower_;
            }

            //- Return coupled interface coefficients, lower
            const FieldField<CoeffField, Type>& coupleLower() const
            {
                return coupleLower_;
            }
			
			//- Access to coupled interfaces
			typename BlockLduInterfaceFieldPtrsList<Type>::Type& interfaces()
			{
				return interfaces_;
			}


        // Matrix structure

            //- Return true if there is a diagonal
            bool thereIsDiag() const
            {
                return (diagPtr_);
            }

            //- Return true if upper triangle is allocated
            bool thereIsUpper() const
            {
                return (upperPtr_);
            }

            //- Return true if lower triangle is allocated
            bool thereIsLower() const
            {
                return (lowerPtr_);
            }

            //- Return true if matrix is empty
            bool empty() const;

            //- Return true if matrix is diagonal-only
            bool diagonal() const;

            //- Return true if matrix is symmetric
            bool symmetric() const;

            //- Return true if matrix is asymmetric
            bool asymmetric() const;

            //- Return true if matrix is component-coupled
            bool componentCoupled() const;


        // Operations

            //- Sum off-diagonal coefficients and add to diagonal
            void sumDiag();

            //- Sum negative off-diagonal coefficients and add to diagonal
            void negSumDiag();

            //- Check matrix for diagonal dominance
            void check() const;

            //- Relax matrix
            void relax
            (
                const TypeField& x,
                TypeField& b,
                const scalar alpha
            );

            //- Matrix multiplication
            void Amul
            (
                TypeField& Ax,
                const TypeField& x
            ) const;

            //- Matrix multiplication without coupled interface update
            void AmulCore
            (
                TypeField& Ax,
                const TypeField& x
            ) const;

            //- Matrix transpose multiplication
            void Tmul
            (
                TypeField& Ax,
                const TypeField& x
            ) const;

            //- Matrix transpose multiplication without coupled
           //  interface update
            void TmulCore
            (
                TypeField& Ax,
                const TypeField& x
            ) const;


            //- Return decoupled b
            void segregateB
            (
                TypeField& sMul,
                const TypeField& x
            ) const;


        // Coupled interface functionality

            //- Initialise the update of coupled interfaces
            //  for Amul operations
            void initInterfaces
            (
                const FieldField<CoeffField, Type>& interfaceCoeffs,
                TypeField& result,
                const TypeField& psi
            ) const;

            //- Update coupled interfaces
            void updateInterfaces
            (
                const FieldField<CoeffField, Type>& interfaceCoeffs,
                TypeField& result,
                const TypeField& psi
            ) const;


        // Constraint manipulation

            //- Set constrained value in a prescribed point
            void setValue
            (
                const label eqnIndex,
                const Type& value
            );


        // Residual calculation

            //- Calculate residual
            tmp<TypeField> residual
            (
                const TypeField& x
            ) const;

            tmp<TypeField> residual
            (
                const TypeField& x,
                const TypeField& b
            ) const;


        // H-operations

            //- Return L-U vector-matrix multiplication in row-form
            tmp<TypeField> H(const TypeField&) const;

            //- Return L-U  vector-matrix multiplication in off-diagonal form
            tmp<TypeField> faceH(const TypeField&) const;


    // Member operators

        void operator=(const BlockLduMatrix<Type>&);

        void negate();

        void operator+=(const BlockLduMatrix<Type>&);
        void operator-=(const BlockLduMatrix<Type>&);

        void operator*=(const scalarField&);
        void operator*=(const scalar);


    // Ostream operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const BlockLduMatrix<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockLduMatrix.C"
#   include "BlockLduMatrixOperations.C"
#   include "BlockLduMatrixUpdateInterfaces.C"
#   include "BlockLduMatrixATmul.C"
#   include "BlockLduMatrixHOps.C"
#   include "BlockLduMatrixDecouple.C"
#   include "BlockLduMatrixDecoupledHOps.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
