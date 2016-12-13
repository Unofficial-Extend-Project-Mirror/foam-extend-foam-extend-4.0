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
    ILUCp

Description
    ILU preconditioning with arbitrary level of fill in (p), based on Crout
    algorithm.

    Reference: Saad, Y.: Iterative Methods for Sparse Linear Systems (2nd
    Edition), SIAM, 2003.

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved

SourceFiles
    ILUCp.C

\*---------------------------------------------------------------------------*/

#ifndef ILUCp_H
#define ILUCp_H

#include "lduMatrix.H"
#include "extendedLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ILUCp Declaration
\*---------------------------------------------------------------------------*/

class ILUCp
:
    public lduPreconditioner
{
    // Private Data

        //- Extended ldu matrix
        extendedLduMatrix extMatrix_;

        //- Preconditioned diagonal
        scalarField preconDiag_;

        //- Temporary working diagonal
        scalar zDiag_;

        //- Temporary working row field
        scalarField z_;

        //- Temporary Working column field
        scalarField w_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ILUCp(const ILUCp&);

        //- Disallow default bitwise assignment
        void operator=(const ILUCp&);

        //- Calculate LU factorization
        void calcFactorization();


public:

    //- Runtime type information
    TypeName("ILUCp");


    // Constructors

        //- Construct from matrix and dictionary
        ILUCp
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~ILUCp();


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt = 0
        ) const;

        //- Execute preconditioning with matrix transpose
        virtual void preconditionT
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
