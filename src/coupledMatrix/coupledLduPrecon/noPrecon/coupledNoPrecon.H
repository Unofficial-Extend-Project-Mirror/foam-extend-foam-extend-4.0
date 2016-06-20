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
    coupledNoPrecon

Description
    No preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    coupledNoPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef coupledNoPrecon_H
#define coupledNoPrecon_H

#include "coupledLduPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class coupledNoPrecon Declaration
\*---------------------------------------------------------------------------*/

class coupledNoPrecon
:
    public coupledLduPrecon
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        coupledNoPrecon(const coupledNoPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const coupledNoPrecon&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        coupledNoPrecon
        (
            const coupledLduMatrix& matrix,
            const PtrList<FieldField<Field, scalar> >& bouCoeffs,
            const PtrList<FieldField<Field, scalar> >& intCoeffs,
            const lduInterfaceFieldPtrsListList& interfaces,
            const dictionary& dict
        )
        :
            coupledLduPrecon
            (
                matrix,
                bouCoeffs,
                intCoeffs,
                interfaces
            )
        {}


    // Destructor

        virtual ~coupledNoPrecon()
        {}


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            FieldField<Field, scalar>& x,
            const FieldField<Field, scalar>& b,
            const direction
        ) const
        {
            x = b;
        }

        //- Execute preconditioning for the transpose matrix
        virtual void preconditionT
        (
            FieldField<Field, scalar>& x,
            const FieldField<Field, scalar>& b,
            const direction
        ) const
        {
            x = b;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
