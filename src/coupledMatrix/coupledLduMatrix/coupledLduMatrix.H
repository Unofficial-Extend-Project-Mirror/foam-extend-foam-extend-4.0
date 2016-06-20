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
    coupledLduMatrix

Description
    Collection of lduMatrices solved together as a block system

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    coupledLduMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef coupledLduMatrix_H
#define coupledLduMatrix_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class coupledLduMatrix Declaration
\*---------------------------------------------------------------------------*/

class coupledLduMatrix
:
    public PtrList<lduMatrix>
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        coupledLduMatrix(const coupledLduMatrix&);

        //- Disallow default bitwise assignment
        void operator=(const coupledLduMatrix&);


public:

    // Static data members

        // Declare name of the class and it's debug switch
        ClassName("coupledLduMatrix");


    // Constructors

        //- Construct given size
        coupledLduMatrix(const label size);


    // Destructor

        ~coupledLduMatrix();


    // Member Functions

        // Matrix structure

            //- Is the matrix diagonal-only?
            bool diagonal() const;

            //- Is the matrix symmetric?
            bool symmetric() const;

            //- Is the matrix asymmetric?
            bool asymmetric() const;


        // Matrix-vector multiplication

            //- Matrix multiplication with updated coupled interfaces.
            void Amul
            (
                FieldField<Field, scalar>& result,
                const FieldField<Field, scalar>& x,
                const PtrList<FieldField<Field, scalar> >& bouCoeffs,
                const lduInterfaceFieldPtrsListList& interfaces,
                const direction cmpt
            ) const;

            //- Matrix transpose multiplication with updated coupled interfaces
            void Tmul
            (
                FieldField<Field, scalar>& result,
                const FieldField<Field, scalar>& x,
                const PtrList<FieldField<Field, scalar> >& intCoeffs,
                const lduInterfaceFieldPtrsListList& interfaces,
                const direction cmpt
            ) const;


            //- Initialise the update of coupled interfaces
            //  for matrix operations
            void initMatrixInterfaces
            (
                const PtrList<FieldField<Field, scalar> >& coupleCoeffs,
                const lduInterfaceFieldPtrsListList& interfaces,
                const FieldField<Field, scalar>& x,
                FieldField<Field, scalar>& result,
                const direction cmpt,
                const bool switchToLhs = false
            ) const;

            //- Update coupled interfaces for matrix operations
            void updateMatrixInterfaces
            (
                const PtrList<FieldField<Field, scalar> >& coupleCoeffs,
                const lduInterfaceFieldPtrsListList& interfaces,
                const FieldField<Field, scalar>& x,
                FieldField<Field, scalar>& result,
                const direction cmpt,
                const bool switchToLhs = false
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
