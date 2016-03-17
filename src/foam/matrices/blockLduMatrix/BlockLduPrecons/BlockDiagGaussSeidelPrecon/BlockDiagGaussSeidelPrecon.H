/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    BlockDiagGaussSeidelPrecon

Description
    Diagonally-corrected block Gauss-Seidel preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    BlockDiagGaussSeidelPrecon.C
    BlockDiagGaussSeidelPrecon.C
    BlockDiagGaussSeidelPreconDecoupled.C

\*---------------------------------------------------------------------------*/

#ifndef BlockDiagGaussSeidelPrecon_H
#define BlockDiagGaussSeidelPrecon_H

#include "BlockLduPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class BlockDiagGaussSeidelPrecon Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockDiagGaussSeidelPrecon
:
    public BlockLduPrecon<Type>
{
    // Private Data

        //- Inverse of diagonal diagonal
        CoeffField<Type> invDiag_;

        //- Off-diag part of diagonal
        CoeffField<Type> LUDiag_;

        //- Temporary space for updated decoupled source
        //  Initialised with zero size and resized on first use
        mutable Field<Type> bPlusLU_;

        //- Temporary space for solution intermediate
        mutable Field<Type> bPrime_;

        //- Number of sweeps
        const label nSweeps_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockDiagGaussSeidelPrecon(const BlockDiagGaussSeidelPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const BlockDiagGaussSeidelPrecon&);


        //- Calculate inverse diagonal
        void calcInvDiag();

        //- Block Gauss-Seidel sweep, symmetric matrix
        template<class DiagType, class ULType>
        void BlockSweep
        (
            Field<Type>& x,
            const Field<DiagType>& dD,
            const Field<ULType>& upper,
            const Field<Type>& b
        ) const;

        //- Block Gauss-Seidel sweep, asymmetric matrix
        template<class DiagType, class ULType>
        void BlockSweep
        (
            Field<Type>& x,
            const Field<DiagType>& dD,
            const Field<ULType>& lower,
            const Field<ULType>& upper,
            const Field<Type>& b
        ) const;


        // Decoupled operations, used in template specialisation

            //- Calculate inverse diagonal, decoupled version
            void calcDecoupledInvDiag();

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
    TypeName("DiagGaussSeidel");


    // Constructors

        //- Construct from matrix for smoother use
        BlockDiagGaussSeidelPrecon
        (
            const BlockLduMatrix<Type>& matrix
        );

        //- Construct from components
        BlockDiagGaussSeidelPrecon
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


    // Destructor

        virtual ~BlockDiagGaussSeidelPrecon()
        {}


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockDiagGaussSeidelPrecon.C"
#   include "BlockDiagGaussSeidelPreconDecoupled.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
