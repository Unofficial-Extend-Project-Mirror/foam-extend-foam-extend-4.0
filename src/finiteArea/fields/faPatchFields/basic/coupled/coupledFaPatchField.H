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
    coupledFaPatchField

Description

SourceFiles
    coupledFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFaPatchField_H
#define coupledFaPatchField_H

#include "lduInterfaceField.H"
#include "faPatchField.H"
#include "coupledFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class coupledFaPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class coupledFaPatchField
:
    public lduInterfaceField,
    public faPatchField<Type>
{

public:

    //- Runtime type information
    TypeName(coupledFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        coupledFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch and internal field and patch field
        coupledFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        coupledFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given coupledFaPatchField onto a new patch
        coupledFaPatchField
        (
            const coupledFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        coupledFaPatchField
        (
            const coupledFaPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<Type> > clone() const = 0;

        //- Construct as copy setting internal field reference
        coupledFaPatchField
        (
            const coupledFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<Type> > clone
        (
            const DimensionedField<Type, areaMesh>&
        ) const = 0;


    // Member functions

        // Access

            //- Return true if this patch field is derived from
            //  coupledFaPatchField<Type>.
            virtual bool coupled() const
            {
                return true;
            }

            //- Return neighbour field of internal field
            virtual tmp<Field<Type> > patchNeighbourField() const = 0;


        // Evaluation functions

            //- Return patch-normal gradient
            virtual tmp<Field<Type> > snGrad() const;

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType
            );

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientBoundaryCoeffs() const;


        // Coupled interface functionality

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction,
                const Pstream::commsTypes commsType
            ) const = 0;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "coupledFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
