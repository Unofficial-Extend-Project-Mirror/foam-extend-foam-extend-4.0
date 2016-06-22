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
    iluSmoother

Description
    ILU smoother

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    iluSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef iluSmoother_H
#define iluSmoother_H

#include "CholeskyPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class iluSmoother Declaration
\*---------------------------------------------------------------------------*/

class iluSmoother
:
    public lduSmoother
{
    // Private Data

        //- Cholesky preconditioner
        CholeskyPrecon precon_;

        //- Correction array
        mutable scalarField xCorr_;

        //- Residual array
        mutable scalarField residual_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        iluSmoother(const iluSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const iluSmoother&);


public:

    //- Runtime type information
    TypeName("ILU");


    // Constructors

        //- Construct from components
        iluSmoother
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );


    // Destructor

        virtual ~iluSmoother()
        {}


    // Member Functions

        //- Execute smoothing
        virtual void smooth
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt,
            const label nSweeps
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
