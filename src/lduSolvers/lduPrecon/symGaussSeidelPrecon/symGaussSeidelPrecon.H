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
    symGaussSeidelPrecon

Description
    Symmetric Gauss-Seidel preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    symGaussSeidelPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef symGaussSeidelPrecon_H
#define symGaussSeidelPrecon_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class symGaussSeidelPrecon Declaration
\*---------------------------------------------------------------------------*/

class symGaussSeidelPrecon
:
    public lduMatrix::preconditioner
{
    // Private Data

        // NOTE: Coupled boundary is treated as an effective Jacobi
        // interface in the boundary.  Note: there is a change of sign
        // in the coupled interface update.  The reason for this is
        // that the internal coefficients are all located at the
        // l.h.s. of the matrix whereas the "implicit" coefficients on
        // the coupled boundaries are all created as if the
        // coefficient contribution is of a source-kind (i.e. they
        // have a sign as if they are on the r.h.s. of the matrix.  To
        // compensate for this, it is necessary to turn the sign of
        // the contribution.

        //- Negative boundary coefficients
        FieldField<Field, scalar> mBouCoeffs_;

        //- Temporary space for solution intermediate
        mutable scalarField bPrime_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        symGaussSeidelPrecon(const symGaussSeidelPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const symGaussSeidelPrecon&);


public:

    //- Runtime type information
    TypeName("SymGaussSeidel");


    // Constructors

        //- Construct from components
        symGaussSeidelPrecon
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );


        //- Construct from dictionary
        symGaussSeidelPrecon
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~symGaussSeidelPrecon()
        {}


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            scalarField& wA,
            const scalarField& rA,
            const direction cmpt = 0
        ) const;

        //- Execute preconditioning on a transpose matrix
        virtual void preconditionT
        (
            scalarField& wT,
            const scalarField& rT,
            const direction cmpt = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
