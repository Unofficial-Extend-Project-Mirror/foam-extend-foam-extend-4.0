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
    amgPrecon

Description
    AMG preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    amgPreconI.H
    amgPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef amgPrecon_H
#define amgPrecon_H

#include "lduMatrix.H"
#include "Switch.H"
#include "amgCycle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class amgPrecon Declaration
\*---------------------------------------------------------------------------*/

class amgPrecon
:
    public lduPreconditioner
{
    // Private Data

        //- Cycle type
        const amgCycle::cycleType cycle_;

        //- Number of pre-sweeps
        const label nPreSweeps_;

        //- Number of post-sweeps
        const label nPostSweeps_;

        //- Maximum number of levels
        const label nMaxLevels_;

        //- Scaling
        const Switch scale_;

        //- AMG cycle pointer
        autoPtr<amgCycle> amgPtr_;

        //- x buffer
        mutable scalarField xBuffer_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        amgPrecon(const amgPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const amgPrecon&);


public:

    //- Runtime type information
    TypeName("AMG");


    // Constructors

        //- Construct from matrix and dictionary
        amgPrecon
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaceFields,
            const dictionary& dict
        );


    // Destructor

        virtual ~amgPrecon();


    // Member Functions

        //- Return number of AMG levels
        label nLevels() const;

        //- Calculate residual
        const scalarField& residual
        (
            const scalarField& x,
            const scalarField& b,
            const direction cmpt
        ) const;

        //- Perform a single cycle
        void cycle
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt
        ) const;

        //- Execute preconditioning
        virtual void precondition
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
