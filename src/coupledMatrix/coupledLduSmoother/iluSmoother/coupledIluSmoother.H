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
    coupledIluSmoother

Description
    ILU smoother for coupled matrix

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    coupledIluSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef coupledIluSmoother_H
#define coupledIluSmoother_H

#include "coupledLduSmoother.H"
#include "coupledCholeskyPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class coupledIluSmoother Declaration
\*---------------------------------------------------------------------------*/

class coupledIluSmoother
:
    public coupledLduSmoother
{
    // Private typedefs

    typedef FieldField<Field, scalar> scalarFieldField;


    // Private Data

        // Ilu preconditioner
        coupledCholeskyPrecon precon_;

        //- Correction array
        mutable scalarFieldField xCorr_;

        //- Residual array
        mutable scalarFieldField residual_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        coupledIluSmoother(const coupledIluSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const coupledIluSmoother&);


public:

    //- Runtime type information
    TypeName("ILU");


    // Constructors

        //- Construct from matrix
        coupledIluSmoother
        (
            const coupledLduMatrix& matrix,
            const PtrList<FieldField<Field, scalar> >& bouCoeffs,
            const PtrList<FieldField<Field, scalar> >& intCoeffs,
            const lduInterfaceFieldPtrsListList& interfaces
        );


    // Destructor

        virtual ~coupledIluSmoother()
        {}


    // Member Functions

        //- Execute smoothing.
        virtual void smooth
        (
            FieldField<Field, scalar>& x,
            const FieldField<Field, scalar>& b,
            const direction cmpt,
            const label nSweeps
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
