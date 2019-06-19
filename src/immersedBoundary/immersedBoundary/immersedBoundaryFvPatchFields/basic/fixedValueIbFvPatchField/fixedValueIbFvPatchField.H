/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::fixedValueIbFvPatchField

Description
    Foam::fixedValueIbFvPatchField

Author
    Hrvoje Jasak

SourceFiles
    fixedValueIbFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValueIbFvPatchField_H
#define fixedValueIbFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "immersedBoundaryFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fixedValueIbFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedValueIbFvPatchField
:
    public fixedValueFvPatchField<Type>,
    public immersedBoundaryFieldBase<Type>
{
    // Private data

        // Fixed value field
        // Note: defining field carry values on faces of the IB patch
        // represented as a triangulated surface

            //- Defining value field
            Field<Type> triValue_;


    // Private Member Functions

        //- Update IB value and gradient
        void updateIbValues();


public:

    //- Runtime type information
    TypeName("fixedValueIb");


    // Constructors

        //- Construct from patch and internal field
        fixedValueIbFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedValueIbFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueIbFvPatchField
        //  onto a new patch
        fixedValueIbFvPatchField
        (
            const fixedValueIbFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedValueIbFvPatchField
        (
            const fixedValueIbFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new fixedValueIbFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedValueIbFvPatchField
        (
            const fixedValueIbFvPatchField<Type>&,
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
                new fixedValueIbFvPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~fixedValueIbFvPatchField()
    {}


    // Member functions

        // Return defining fields
        // Note: defining fields carry values on faces of the IB patch
        // represented as a triangulated surface

            //- Return reference value
            const Field<Type>& triValue() const
            {
                return triValue_;
            }

            //- Return access to reference value
            Field<Type>& triValue()
            {
                return triValue_;
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

            //- Update on mesh motion
            virtual void updateOnMotion();


        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );

            //- Manipulate a matrix
            virtual void manipulateMatrix(fvMatrix<Type>& matrix);


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fixedValueIbFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
