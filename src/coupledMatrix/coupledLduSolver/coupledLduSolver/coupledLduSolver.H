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
    coupledLduSolver

Description
    Coupled ldu solver virtual base class

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    coupledLduSolver.C

\*---------------------------------------------------------------------------*/

#ifndef coupledLduSolver_H
#define coupledLduSolver_H

#include "coupledLduMatrix.H"
#include "coupledSolverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class coupledLduSolver Declaration
\*---------------------------------------------------------------------------*/

class coupledLduSolver
{
    // Private data

        //- Field name
        word fieldName_;


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

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            coupledLduSolver,
            symMatrix,
            (
                const word& fieldName,
                const coupledLduMatrix& matrix,
                const PtrList<FieldField<Field, scalar> >& bouCoeffs,
                const PtrList<FieldField<Field, scalar> >& intCoeffs,
                const lduInterfaceFieldPtrsListList& interfaces,
                const dictionary& solverData
            ),
            (
                fieldName,
                matrix,
                bouCoeffs,
                intCoeffs,
                interfaces,
                solverData
            )
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            coupledLduSolver,
            asymMatrix,
            (
                const word& fieldName,
                const coupledLduMatrix& matrix,
                const PtrList<FieldField<Field, scalar> >& bouCoeffs,
                const PtrList<FieldField<Field, scalar> >& intCoeffs,
                const lduInterfaceFieldPtrsListList& interfaces,
                const dictionary& solverData
            ),
            (
                fieldName,
                matrix,
                bouCoeffs,
                intCoeffs,
                interfaces,
                solverData
            )
        );


        // Constructors

            coupledLduSolver
            (
                const word& fieldName,
                const coupledLduMatrix& matrix,
                const PtrList<FieldField<Field, scalar> >& bouCoeffs,
                const PtrList<FieldField<Field, scalar> >& intCoeffs,
                const lduInterfaceFieldPtrsListList& interfaces
            )
            :
                fieldName_(fieldName),
                matrix_(matrix),
                bouCoeffs_(bouCoeffs),
                intCoeffs_(intCoeffs),
                interfaces_(interfaces)
            {}


        // Selectors

            //- Return a new solver
            static autoPtr<coupledLduSolver> New
            (
                const word& fieldName,
                const coupledLduMatrix& matrix,
                const PtrList<FieldField<Field, scalar> >& bouCoeffs,
                const PtrList<FieldField<Field, scalar> >& intCoeffs,
                const lduInterfaceFieldPtrsListList& interfaces,
                const dictionary& solverData
            );


        // Destructor

            virtual ~coupledLduSolver()
            {}


        // Member functions

            //- Return field name
            const word& fieldName() const
            {
                return fieldName_;
            }

            //- Solve
            virtual coupledSolverPerformance solve
            (
                FieldField<Field, scalar>& x,
                const FieldField<Field, scalar>& b,
                const direction cmpt = 0
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
