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
    Foam::fixedInternalValueFvPatchField

Description
    Boundary condition providing mechanism to set boundary (cell) values
    directly into a matrix, i.e. to set a constraint condition. Default
    behaviour is to act as a zero gradient condition.

SourceFiles
    fixedInternalValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedInternalValueFvPatchField_H
#define fixedInternalValueFvPatchField_H

#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class fixedInternalValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedInternalValueFvPatchField
:
    public zeroGradientFvPatchField<Type>
{
    // Private data

        //- Value to be set in internal field
        Field<Type> refValue_;


    // Private Member Functions

        //- Set refValue in internal field
        void setInInternalField() const;


public:

    //- Runtime type information
    TypeName("fixedInternalValue");


    // Constructors

        //- Construct from patch and internal field
        fixedInternalValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedInternalValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given fixedInternalValueFvPatchField<Type>
        //  onto a new patch
        fixedInternalValueFvPatchField
        (
            const fixedInternalValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedInternalValueFvPatchField
        (
            const fixedInternalValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new fixedInternalValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedInternalValueFvPatchField
        (
            const fixedInternalValueFvPatchField<Type>&,
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
                new fixedInternalValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return reference field value
            const Field<Type>& refValue() const
            {
                return refValue_;
            }

            //- Return access to reference field value
            Field<Type>& refValue()
            {
                return refValue_;
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


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field, sets Updated to false
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );

            //- Manipulate a matrix
            virtual void manipulateMatrix(fvMatrix<Type>& matrix);


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fixedInternalValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
