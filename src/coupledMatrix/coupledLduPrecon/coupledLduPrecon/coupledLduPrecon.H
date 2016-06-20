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
    coupledLduPrecon

Description
    Coupled LDU matrix preconditioner virtual base class

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    coupledLduPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef coupledLduPrecon_H
#define coupledLduPrecon_H

#include "coupledLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class coupledLduPrecon Declaration
\*---------------------------------------------------------------------------*/

class coupledLduPrecon
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        coupledLduPrecon(const coupledLduPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const coupledLduPrecon&);


protected:

    // Protected data

        //- Matrix reference
        const coupledLduMatrix& matrix_;

        //- Boundary coefficients
        const PtrList<FieldField<Field, scalar> >& bouCoeffs_;

        //- Internal coefficients
        const PtrList<FieldField<Field, scalar> >& intCoeffs_;

        //- Interfaces
        const lduInterfaceFieldPtrsListList& interfaces_;


public:

    //- Find the preconditioner name (directly or from a sub-dictionary)
    static word getName(const dictionary&);

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            coupledLduPrecon,
            dictionary,
            (
                const coupledLduMatrix& matrix,
                const PtrList<FieldField<Field, scalar> >& bouCoeffs,
                const PtrList<FieldField<Field, scalar> >& intCoeffs,
                const lduInterfaceFieldPtrsListList& interfaces,
                const dictionary& dict
            ),
            (
                matrix,
                bouCoeffs,
                intCoeffs,
                interfaces,
                dict
            )
        );


    // Constructors

        //- Construct from components
        coupledLduPrecon
        (
            const coupledLduMatrix& matrix,
            const PtrList<FieldField<Field, scalar> >& bouCoeffs,
            const PtrList<FieldField<Field, scalar> >& intCoeffs,
            const lduInterfaceFieldPtrsListList& interfaces
        )
        :
            matrix_(matrix),
            bouCoeffs_(bouCoeffs),
            intCoeffs_(intCoeffs),
            interfaces_(interfaces)
        {}


    // Selectors

        //- Select given matrix data and dictionary
        static autoPtr<coupledLduPrecon> New
        (
            const coupledLduMatrix& matrix,
            const PtrList<FieldField<Field, scalar> >& bouCoeffs,
            const PtrList<FieldField<Field, scalar> >& intCoeffs,
            const lduInterfaceFieldPtrsListList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~coupledLduPrecon()
        {}


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            FieldField<Field, scalar>& x,
            const FieldField<Field, scalar>& b,
            const direction cmpt = 0
        ) const = 0;

        //- Execute preconditioning for the transpose matrix
        virtual void preconditionT
        (
            FieldField<Field, scalar>& x,
            const FieldField<Field, scalar>& b,
            const direction cmpt = 0
        ) const
        {
            notImplemented
            (
                type() +"::preconditionT"
                "(FieldField<Field, scalar>& x, "
                "const FieldField<Field, scalar>& b,"
                "const direction cmpt) const"
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
