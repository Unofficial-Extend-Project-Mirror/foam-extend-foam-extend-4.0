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
    CholeskyPrecon

Description
    Incmplete Cholesky preconditioning with no fill-in for symmeytric matrices

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    CholeskyPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef CholeskyPrecon_H
#define CholeskyPrecon_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CholeskyPrecon Declaration
\*---------------------------------------------------------------------------*/

class CholeskyPrecon
:
    public lduPreconditioner
{
    // Private Data

        //- Preconditioned diagonal
        scalarField preconDiag_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CholeskyPrecon(const CholeskyPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const CholeskyPrecon&);

        //- Precondition the diagonal
        void calcPreconDiag();


public:

    //- Runtime type information
    TypeName("Cholesky");


    // Constructors

        //- Construct from matrix and dictionary
        CholeskyPrecon
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );

        //- Construct from matrix as a smoother
        CholeskyPrecon
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );


    // Destructor

        virtual ~CholeskyPrecon();


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            scalarField& wA,
            const scalarField& rA,
            const direction cmpt = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
