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
    coupledLduSmoother

Description
    Coupled LDU matrix smoother virtual base class.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    coupledLduSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef coupledLduSmoother_H
#define coupledLduSmoother_H

#include "coupledLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class coupledLduSmoother Declaration
\*---------------------------------------------------------------------------*/

class coupledLduSmoother
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        coupledLduSmoother(const coupledLduSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const coupledLduSmoother&);


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

    //- Find the smoother name (directly or from a sub-dictionary)
    static word getName(const dictionary&);

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            coupledLduSmoother,
            word,
            (
                const coupledLduMatrix& matrix,
                const PtrList<FieldField<Field, scalar> >& bouCoeffs,
                const PtrList<FieldField<Field, scalar> >& intCoeffs,
                const lduInterfaceFieldPtrsListList& interfaces
            ),
            (
                matrix,
                bouCoeffs,
                intCoeffs,
                interfaces
            )
        );


    // Constructors

        //- Construct from components
        coupledLduSmoother
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

        //- Select given matrix data and type
        static autoPtr<coupledLduSmoother> New
        (
            const coupledLduMatrix& matrix,
            const PtrList<FieldField<Field, scalar> >& bouCoeffs,
            const PtrList<FieldField<Field, scalar> >& intCoeffs,
            const lduInterfaceFieldPtrsListList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~coupledLduSmoother()
        {}


    // Member Functions

        //- Execute smoothing.
        virtual void smooth
        (
            FieldField<Field, scalar>& x,
            const FieldField<Field, scalar>& b,
            const direction cmpt,
            const label nSweeps
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
