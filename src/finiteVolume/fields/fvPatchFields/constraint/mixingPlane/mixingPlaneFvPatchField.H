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
    mixingPlaneFvPatchField

Description
    MixingPlane interface patch field.

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    mixingPlaneFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPlaneFvPatchField_H
#define mixingPlaneFvPatchField_H

#include "coupledFvPatchField.H"
#include "mixingPlaneLduInterfaceField.H"
#include "mixingPlaneFvPatch.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class mixingPlaneFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mixingPlaneFvPatchField
:
    public mixingPlaneLduInterfaceField,
    public coupledFvPatchField<Type>
{
    // Private data

        //- Local reference cast into the mixingPlane patch
        const mixingPlaneFvPatch& mixingPlanePatch_;

        //- Interpolation type
        mutable mixingPlaneInterpolation::mixingType mixing_;

        //- Name of flux field
        word phiName_;

        //- Flux mask for averaging.  Size is equal to the number of bands
        //  Flux mask behaves like valueFraction:
        //      1 = use interpolated value from mixing plane interpolator
        //      0 = use patch internal field to mimic zero gradient behaviour
        mutable scalarField fluxMask_;

        //- Flux weights for averaging
        mutable scalarField fluxWeights_;


    // Private member functions

        //- Read mixing type
        void readMixingType() const;

        //- Calculate flux mask and weights
        void calcFluxMask() const;


public:

    //- Runtime type information
    TypeName(mixingPlaneFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        mixingPlaneFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mixingPlaneFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixingPlaneFvPatchField onto a new patch
        mixingPlaneFvPatchField
        (
            const mixingPlaneFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mixingPlaneFvPatchField
        (
            const mixingPlaneFvPatchField&
        );

        //- Construct as copy setting internal field reference
        mixingPlaneFvPatchField
        (
            const mixingPlaneFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new mixingPlaneFvPatchField<Type>(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new mixingPlaneFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

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


        // Access

            //- Return reference to mixingPlanePatch patch
            const mixingPlaneFvPatch& mixingPlanePatch() const
            {
                return mixingPlanePatch_;
            }

            //- Return shadow patch field
            const mixingPlaneFvPatchField<Type>& shadowPatchField() const;

            //- Return reference to mixing patch
            const mixingPlaneInterpolation::mixingType& mixing() const;

            //- Return flux mask.  Flux mask is calculated from master fluxes
            //  and relates to the mixing plane profile
            const scalarField& fluxMask() const;

            //- Return flux weights.  Flux weights relate to local faces
            //  and are used in weighted interpolation
            const scalarField& fluxWeights() const;


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return neighbour patch field
            virtual tmp<Field<Type> > patchNeighbourField() const;

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType
            );


        // Mixing plane discretisation support

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

//             //- Return the matrix diagonal coefficients corresponding to the
//             //  evaluation of the gradient of this patchField
//             virtual tmp<Field<Type> > gradientInternalCoeffs() const;

//             //- Return the matrix source coefficients corresponding to the
//             //  evaluation of the gradient of this patchField
//             virtual tmp<Field<Type> > gradientBoundaryCoeffs() const;


        // Patch interpolation and patch flux

            //- Calculate patch face interpolate given weighting factors
            virtual void patchInterpolate
            (
                GeometricField<Type, fvsPatchField, surfaceMesh>& fField,
                const scalarField& pL
            ) const;

            //- Calculate patch face interpolate given two weighting factors
            virtual void patchInterpolate
            (
                GeometricField<Type, fvsPatchField, surfaceMesh>& fField,
                const scalarField& pL,
                const scalarField& pY
            ) const;

            //- Calculate patch flux
            virtual void patchFlux
            (
                GeometricField<Type, fvsPatchField, surfaceMesh>& flux,
                const fvMatrix<Type>& matrix
            ) const;


        // Matrix manipulation

            //- Transform given patch component field
            virtual void transformCoupleField
            (
                scalarField& f,
                const direction cmpt
            ) const
            {
                mixingPlaneLduInterfaceField::transformCoupleField(f,cmpt);
            }

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;


        // MixingPlane coupled interface functions

            //- Does the patch field perform the transformation
            virtual bool doTransform() const
            {
                return
               !(
                    mixingPlanePatch_.parallel()
                 || pTraits<Type>::rank == 0
                );
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return mixingPlanePatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return mixingPlanePatch_.reverseT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<Type>::rank;
            }

             //- Write
             virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mixingPlaneFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
