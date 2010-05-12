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
    fixedGradientFaPatchField

Description

SourceFiles
    fixedGradientFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedGradientFaPatchField_H
#define fixedGradientFaPatchField_H

#include "faPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class fixedGradientFaPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedGradientFaPatchField
:
    public faPatchField<Type>
{
    // Private data

        Field<Type> gradient_;


public:

    //- Runtime type information
    TypeName("fixedGradient");


    // Constructors

        //- Construct from patch and internal field
        fixedGradientFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedGradientFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given fixedGradientFaPatchField<Type>
        //  onto a new patch
        fixedGradientFaPatchField
        (
            const fixedGradientFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        fixedGradientFaPatchField
        (
            const fixedGradientFaPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<Type> > clone() const
        {
            return tmp<faPatchField<Type> >
            (
                new fixedGradientFaPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedGradientFaPatchField
        (
            const fixedGradientFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<Type> > clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return tmp<faPatchField<Type> >
            (
                new fixedGradientFaPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Return defining fields

            //- Return gradient at boundary
            virtual Field<Type>& gradient()
            {
                return gradient_;
            }

            virtual const Field<Type>& gradient() const
            {
                return gradient_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const faPatchFieldMapper&
            );

            //- Reverse map the given faPatchField onto this faPatchField
            virtual void rmap
            (
                const faPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type> > snGrad() const
            {
                return gradient_;
            }

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
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


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fixedGradientFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
