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
    Foam::FDICPreconditioner

Description
    Faster version of the DICPreconditioner diagonal-based incomplete
    Cholesky preconditioner for symmetric matrices
    (symmetric equivalent of DILU) in which the the reciprocal of the
    preconditioned diagonal and the upper coefficients divided by the diagonal
    are calculated and stored.

SourceFiles
    FDICPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef FDICPreconditioner_H
#define FDICPreconditioner_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class FDICPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class FDICPreconditioner
:
    public lduPreconditioner
{
    // Private data

        //- The reciprocal preconditioned diagonal
        scalarField rD_;
        scalarField rDuUpper_;
        scalarField rDlUpper_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        FDICPreconditioner(const FDICPreconditioner&);

        //- Disallow default bitwise assignment
        void operator=(const FDICPreconditioner&);


public:

    //- Runtime type information
    TypeName("FDIC");


    // Constructors

        //- Construct from matrix components and preconditioner solver controls
        FDICPreconditioner
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~FDICPreconditioner()
        {}


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
