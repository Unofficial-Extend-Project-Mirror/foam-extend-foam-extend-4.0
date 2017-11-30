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
    Foam::immersedBoundaryFvPatchField

Description
    Foam::immersedBoundaryFvPatchField

Author
    Hrvoje Jasak

SourceFiles
    immersedBoundaryFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryFvPatchField_H
#define immersedBoundaryFvPatchField_H

#include "fvPatchField.H"
#include "immersedBoundaryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class immersedBoundaryFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class immersedBoundaryFvPatchField
:
    public fvPatchField<Type>
{
    // Private data

        //- Local reference to immersed boundary patch
        const immersedBoundaryFvPatch& ibPatch_;

        //- Local reference to fvMesh
        const fvMesh& mesh_;

        // Defining fields
        // Note: defining fields carry values on faces of the IB patch
        // represented as a triangulated surface

            //- Defining value field
            Field<Type> refValue_;

            //- Defining normal gradient field
            Field<Type> refGrad_;

            //- Fraction (0-1) of value used for boundary condition
            scalarField refValueFraction_;


        // Dead cell controls

            //- Set dead cell value
            Switch setDeadCellValue_;

            //- Dead cell value
            Type deadCellValue_;


        // Fields on the IB face centres
        // Field size equals the number of cells intesected with the IB

            //- Time index for last update of mesh or moving boundary
            mutable label ibUpdateTimeIndex_;

            //- Value field
            mutable Field<Type> ibValue_;

            //- Normal gradient field
            mutable Field<Type> ibGrad_;

            //- Value fraction
            mutable scalarField ibValueFraction_;


    // Private Member Functions

        //- Update IB value and gradient
        void updateIbValues() const;

        //- Impose condition in dead cells
        void imposeDeadCondition();


        // Manipulate matrix

            //- Check and correct zero diagonal in dead cells
            void correctDiag
            (
                fvMatrix<Type>& eqn
            ) const;

            //- Impose fixed gradient condition by manipulating matrix
            //  Note: reconsider pre-factor for diffusivity
            //  HJ, 16/Apr/2012
            void correctOffDiag
            (
                fvMatrix<Type>& eqn
            ) const;


protected:

    // Protected Member Functions

        //- Does immersed boundary patch field need motion update?
        bool motionUpdateRequired() const;

        //- Execute immersed boundary patch field motion update
        virtual void motionUpdate() const;


public:

    //- Runtime type information
    TypeName(immersedBoundaryFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        immersedBoundaryFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedBoundaryFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given immersedBoundaryFvPatchField
        //  onto a new patch
        immersedBoundaryFvPatchField
        (
            const immersedBoundaryFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        immersedBoundaryFvPatchField
        (
            const immersedBoundaryFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new immersedBoundaryFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        immersedBoundaryFvPatchField
        (
            const immersedBoundaryFvPatchField<Type>&,
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
                new immersedBoundaryFvPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~immersedBoundaryFvPatchField()
    {}


    // Member functions

        //- Return reference to immersed boundary patch
        const immersedBoundaryFvPatch& ibPatch() const
        {
            return ibPatch_;
        }


        // Return defining fields
        // Note: defining fields carry values on faces of the IB patch
        // represented as a triangulated surface

            //- Return reference value
            const Field<Type>& refValue() const
            {
                return refValue_;
            }

            //- Return access to reference value
            Field<Type>& refValue()
            {
                return refValue_;
            }

            //- Return reference gradient
            const Field<Type>& refGrad() const
            {
                return refGrad_;
            }

            //- Return access to reference gradient
            Field<Type>& refGrad()
            {
                return refGrad_;
            }

            //- Return value fraction
            virtual const scalarField& valueFraction() const
            {
                return fixesValue_;
            }

            //- Return access to value fraction
            scalarField& valueFraction()
            {
                return valueFraction_;
            }


        //- Return fields on intersection points interacting with the IB

            //- Return IB field
            const Field<Type>& ibValue() const;

            //- Return IB surface-normal gradient
            const Field<Type>& ibGrad() const;

            //- Return IB cell field
            tmp<Field<Type> > ibCellValue() const;

            //- Return IB sample point field
            tmp<Field<Type> > ibSamplingPointValue() const;


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            )
            {}

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            )
            {}


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            void updateCoeffs();

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const
            {
                return tmp<Field<Type> >(new Field<Type>(0));
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const
            {
                return tmp<Field<Type> >(new Field<Type>(0));
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            tmp<Field<Type> > gradientInternalCoeffs() const
            {
                return tmp<Field<Type> >(new Field<Type>(0));
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            tmp<Field<Type> > gradientBoundaryCoeffs() const
            {
                return tmp<Field<Type> >(new Field<Type>(0));
            }

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<Type>& matrix);


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "immersedBoundaryFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
