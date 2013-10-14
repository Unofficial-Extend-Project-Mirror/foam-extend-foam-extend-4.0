/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    newDirectionMixedFvPatchField

Description
    Doubly mixed fixed value-fixed gradient boundary condition
    separated into a normal and a tangential component given a
    direction vector.  The mixture is controlled by two separate
    valueFraction coefficients in the normal and tangential direction.

SourceFiles
    newDirectionMixedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef newDirectionMixedFvPatchField_H
#define newDirectionMixedFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class newDirectionMixedFvPatchField
:
    public fvPatchField<Type>
{
    // Private data

        //- Value field
        Field<Type> refValue_;

        //- Gradient field
        Field<Type> refGrad_;

        //- Normal direction
        vectorField nHat_;

        //- Fraction (0-1) of value used for normal component
        scalarField normalValueFraction_;

        //- Fraction (0-1) of value used for tangential component
        scalarField tangentialValueFraction_;


    // Private Member Functions

        //- Check and normalize nHat
        void checkNHat();


public:

    //- Runtime type information
    TypeName("newDirectionMixed");


    // Constructors

        //- Construct from patch and internal field
        newDirectionMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        newDirectionMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given newDirectionMixedFvPatchField onto
        //  a new patch
        newDirectionMixedFvPatchField
        (
            const newDirectionMixedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new newDirectionMixedFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        newDirectionMixedFvPatchField
        (
            const newDirectionMixedFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new newDirectionMixedFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return true;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Return defining fields

            virtual Field<Type>& refValue()
            {
                return refValue_;
            }

            virtual const Field<Type>& refValue() const
            {
                return refValue_;
            }

            virtual Field<Type>& refGrad()
            {
                return refGrad_;
            }

            virtual const Field<Type>& refGrad() const
            {
                return refGrad_;
            }

            virtual vectorField& nHat()
            {
                return nHat_;
            }

            virtual const vectorField& nHat() const
            {
                return nHat_;
            }

            virtual scalarField& normalValueFraction()
            {
                return normalValueFraction_;
            }

            virtual const scalarField& normalValueFraction() const
            {
                return normalValueFraction_;
            }


            virtual scalarField& tangentialValueFraction()
            {
                return tangentialValueFraction_;
            }

            virtual const scalarField& tangentialValueFraction() const
            {
                return tangentialValueFraction_;
            }


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type> > snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate();

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


    // Member operators

        virtual void operator=(const fvPatchField<Type>&) {}
        virtual void operator+=(const fvPatchField<Type>&) {}
        virtual void operator-=(const fvPatchField<Type>&) {}
        virtual void operator*=(const fvPatchField<Type>&) {}
        virtual void operator/=(const fvPatchField<Type>&) {}

        virtual void operator=(const Field<Type>&) {}
        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}
        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "newDirectionMixedFvPatchField.C"
#else
#   ifdef xlC
#   pragma implementation("newDirectionMixedFvPatchField.C")
#   endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
