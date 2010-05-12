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
    coupledCholeskyPrecon

Description
    Incomplete Cholesky preconditioning with no fill-in.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    coupledCholeskyPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef coupledCholeskyPrecon_H
#define coupledCholeskyPrecon_H

#include "coupledLduPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class coupledCholeskyPrecon Declaration
\*---------------------------------------------------------------------------*/

class coupledCholeskyPrecon
:
    public coupledLduPrecon
{
    // Private Data

        //- Preconditioned diagonal
        FieldField<Field, scalar> preconDiag_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        coupledCholeskyPrecon(const coupledCholeskyPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const coupledCholeskyPrecon&);

        //- Precondition the diagonal
        void calcPreconDiag();


public:

    //- Runtime type information
    TypeName("Cholesky");


    // Constructors

        //- Construct from components
        coupledCholeskyPrecon
        (
            const coupledLduMatrix& matrix,
            const PtrList<FieldField<Field, scalar> >& bouCoeffs,
            const PtrList<FieldField<Field, scalar> >& intCoeffs,
            const lduInterfaceFieldPtrsListList& interfaces
        );

        //- Construct from dictionary
        coupledCholeskyPrecon
        (
            const coupledLduMatrix& matrix,
            const PtrList<FieldField<Field, scalar> >& bouCoeffs,
            const PtrList<FieldField<Field, scalar> >& intCoeffs,
            const lduInterfaceFieldPtrsListList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~coupledCholeskyPrecon()
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
