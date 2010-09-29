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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    BlockGaussSeidelPrecon

Description
    Gauss-Seidel preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    BlockGaussSeidelPrecon.C
    BlockGaussSeidelPrecon.C
    BlockGaussSeidelPreconDecoupled.C

\*---------------------------------------------------------------------------*/

#ifndef BlockGaussSeidelPrecon_H
#define BlockGaussSeidelPrecon_H

#include "BlockLduPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BlockGaussSeidelPrecon Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockGaussSeidelPrecon
:
    public BlockLduPrecon<Type>
{
    // Private Data

        //- Temporary space for solution intermediate
        mutable Field<Type> bPrime_;

        //- Number of sweeps
        const label nSweeps_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockGaussSeidelPrecon(const BlockGaussSeidelPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const BlockGaussSeidelPrecon&);


        // Block Gauss-Seidel sweep, symetric matrix
        template<class DiagType, class ULType>
        void BlockSweep
        (
            Field<Type>& x,
            const Field<DiagType>& dD,
            const Field<ULType>& upper,
            const Field<Type>& b
        ) const;

        // Block Gauss-Seidel sweep, asymmetric matrix
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

            //- Execute preconditioning, decoupled version
            void decoupledPrecondition
            (
                Field<Type>& x,
                const Field<Type>& b
            ) const;

            //- Execute preconditioning with matrix transpose, decoupled version
            void decoupledPreconditionT
            (
                Field<Type>& xT,
                const Field<Type>& bT
            ) const;


public:

    //- Runtime type information
    TypeName("GaussSeidel");


    // Constructors

        //- Construct from matrix for smoother use
        BlockGaussSeidelPrecon
        (
            const BlockLduMatrix<Type>& matrix
        )
        :
            BlockLduPrecon<Type>(matrix),
            bPrime_(matrix.lduAddr().size()),
            nSweeps_(1)
        {}

        //- Construct from components
        BlockGaussSeidelPrecon
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        )
        :
            BlockLduPrecon<Type>(matrix),
            bPrime_(matrix.lduAddr().size()),
            nSweeps_(readLabel(dict.lookup("nSweeps")))
        {}


    // Destructor

        virtual ~BlockGaussSeidelPrecon()
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
#   include "BlockGaussSeidelPrecon.C"
#   include "BlockGaussSeidelPreconDecoupled.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
