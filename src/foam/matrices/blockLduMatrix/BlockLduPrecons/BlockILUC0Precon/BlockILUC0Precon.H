/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    BlockILUC0Precon

Description
    Block variant of ILU preconditioning with arbitrary level of fill in (p),
    based on Crout algorithm.

    Reference: Saad, Y.: Iterative Methods for Sparse Linear Systems (2nd
    Edition), SIAM, 2003.

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

SourceFiles
    BlockILUC0Precon.C
    BlockILUC0PreconDecoupled.C

\*---------------------------------------------------------------------------*/

#ifndef BlockILUC0Precon_H
#define BlockILUC0Precon_H

#include "BlockLduPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BlockILUC0Precon Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockILUC0Precon
:
    public BlockLduPrecon<Type>
{
    // Private Data

        //- Preconditioned diagonal
        mutable CoeffField<Type> preconDiag_;

        //- Preconditioned lower
        mutable CoeffField<Type> preconLower_;

        //- Preconditioned upper
        mutable CoeffField<Type> preconUpper_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockILUC0Precon(const BlockILUC0Precon&);

        //- Disallow default bitwise assignment
        void operator=(const BlockILUC0Precon&);

        //- Calculate active type dependant factorization.
        //  Note: both lower, diag and upper have to be of same type
        template<class LDUType>
        void calcActiveTypeFactorization
        (
            Field<LDUType>& preconD,
            Field<LDUType>& preconUpper,
            Field<LDUType>& preconLower
        ) const;

        //- Calculate LU factorization (constructor helper)
        void calcFactorization() const;

        //- Performs forward/backward LU substitution
        //  Note: both lower, diag and upper have to be of same type
        template<class LDUType>
        void LUSubstitute
        (
            Field<Type>& x,
            const Field<LDUType>& preconD,
            const Field<LDUType>& upper,
            const Field<LDUType>& lower,
            const Field<Type>& b
        ) const;

        //- Performs forward/backward LU substitution on transpose system
        //  Note: both lower, diag and upper have to be of same type
        template<class LDUType>
        void LUSubstituteT
        (
            Field<Type>& xT,
            const Field<LDUType>& preconD,
            const Field<LDUType>& upper,
            const Field<LDUType>& lower,
            const Field<Type>& bT
        ) const;


        // Decoupled operations, used in template specialisation

            //- Execute preconditioning, decoupled version
            void decoupledPrecondition
            (
                Field<Type>& x,
                const Field<Type>& b
            ) const;

            //- Execute preconditioning with matrix transpose,
            //  decoupled version
            void decoupledPreconditionT
            (
                Field<Type>& xT,
                const Field<Type>& bT
            ) const;


public:

    //- Runtime type information
    TypeName("ILUC0");


    // Constructors

        //- Construct from components
        BlockILUC0Precon
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


    //- Destructor
    virtual ~BlockILUC0Precon();


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            Field<Type>& x,
            const Field<Type>& b
        ) const;

        //- Execute preconditioning with matrix transpose
        virtual void preconditionT
        (
            Field<Type>& xT,
            const Field<Type>& bT
        ) const;

        //- Re-initialise preconditioner after matrix coefficient update
        virtual void initMatrix();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockILUC0Precon.C"
#   include "BlockILUC0PreconDecoupled.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
