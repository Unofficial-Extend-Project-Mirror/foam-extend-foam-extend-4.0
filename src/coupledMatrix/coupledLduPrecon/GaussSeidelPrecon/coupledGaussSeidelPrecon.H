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
    coupledGaussSeidelPrecon

Description
    Symmetric Gauss-Seidel preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    coupledGaussSeidelPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef coupledGaussSeidelPrecon_H
#define coupledGaussSeidelPrecon_H

#include "coupledLduPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class coupledGaussSeidelPrecon Declaration
\*---------------------------------------------------------------------------*/

class coupledGaussSeidelPrecon
:
    public coupledLduPrecon
{
    // Private Data

        //- Negative boundary coefficients
        PtrList<FieldField<Field, scalar> > mBouCoeffs_;

        //- Temporary space for solution intermediate
        mutable FieldField<Field, scalar> bPrime_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        coupledGaussSeidelPrecon(const coupledGaussSeidelPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const coupledGaussSeidelPrecon&);


        //- Forward sweep on a single matrix
        void forwardSweep
        (
            const lduMatrix& matrix,
            scalarField& x,
            scalarField& bPrime
        ) const;

        //- Reverse sweep on a single matrix
        void reverseSweep
        (
            const lduMatrix& matrix,
            scalarField& x,
            scalarField& bPrime
        ) const;


        //- Forward transpose sweep on a single matrix
        void forwardSweepTranspose
        (
            const lduMatrix& matrix,
            scalarField& x,
            scalarField& bPrime
        ) const;

        //- Reverse transpose sweep on a single matrix
        void reverseSweepTranspose
        (
            const lduMatrix& matrix,
            scalarField& x,
            scalarField& bPrime
        ) const;


public:

    //- Runtime type information
    TypeName("GaussSeidel");


    // Constructors

        //- Construct from components
        coupledGaussSeidelPrecon
        (
            const coupledLduMatrix& matrix,
            const PtrList<FieldField<Field, scalar> >& bouCoeffs,
            const PtrList<FieldField<Field, scalar> >& intCoeffs,
            const lduInterfaceFieldPtrsListList& interfaces
        );

        //- Construct from dictionary
        coupledGaussSeidelPrecon
        (
            const coupledLduMatrix& matrix,
            const PtrList<FieldField<Field, scalar> >& bouCoeffs,
            const PtrList<FieldField<Field, scalar> >& intCoeffs,
            const lduInterfaceFieldPtrsListList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~coupledGaussSeidelPrecon()
        {}


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            FieldField<Field, scalar>& x,
            const FieldField<Field, scalar>& b,
            const direction cmpt = 0
        ) const;

        //- Execute preconditioning for the transpose matrix
        virtual void preconditionT
        (
            FieldField<Field, scalar>& x,
            const FieldField<Field, scalar>& b,
            const direction cmpt = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
