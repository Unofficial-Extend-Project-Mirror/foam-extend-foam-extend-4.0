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
    Foam::noPreconditioner

Description
    Null preconditioner for both symmetric and asymmetric matrices.

SourceFiles
    noPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef noPreconditioner_H
#define noPreconditioner_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class noPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class noPreconditioner
:
    public lduPreconditioner
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        noPreconditioner(const noPreconditioner&);

        //- Disallow default bitwise assignment
        void operator=(const noPreconditioner&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from matrix components and preconditioner solver controls
        noPreconditioner
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~noPreconditioner()
        {}


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            scalarField& wA,
            const scalarField& rA,
            const direction cmpt = 0
        ) const;

        //- Execute preconditioning with matrix transpose
        virtual void preconditionT
        (
            scalarField& wT,
            const scalarField& rT,
            const direction cmpt = 0
        ) const
        {
            precondition(wT, rT, cmpt);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
