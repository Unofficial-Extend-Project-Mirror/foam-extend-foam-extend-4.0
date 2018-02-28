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
    Foam::oversetFvPatchField

Description
    Foam::oversetFvPatchField

SourceFiles
    oversetFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef oversetFvPatchField_H
#define oversetFvPatchField_H

#include "coupledFvPatchField.H"
#include "oversetLduInterfaceField.H"
#include "oversetFvPatch.H"
#include "oversetMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class oversetFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class oversetFvPatchField
:
    public oversetLduInterfaceField,
    public coupledFvPatchField<Type>
{
    // Private data

        //- Local reference cast into the overset patch
        const oversetFvPatch& oversetPatch_;


        // Hole cell controls

            //- Coupled fringe
            Switch coupledFringe_;

            //- Correct interpolation for fringe conservation (use for pressure
            //  only)
            Switch conservativeCorrection_;

            //- Set hole cell value
            Switch setHoleCellValue_;

            //- Hole cell value
            Type holeCellValue_;


        // Fringe matrix coefficients, collected on matrix completion
        // Coefficients must be stored because the matrix only updates
        // live cells, while the acceptor values are updated by interpolation
        // HJ, 21/May/2013

            //- Fringe upper coefficients
            mutable scalarField fringeUpperCoeffs_;

            //- Fringe lower coefficients
            mutable scalarField fringeLowerCoeffs_;


    // Private Member Functions

        //- Set value in hole cells
        template<class Type2>
        void setHoleValues(Field<Type2>& f) const;

        //- Set value in acceptor cells.  Uses interpolation
        template<class Type2>
        void setAcceptorValues(Field<Type2>& f) const;


        // Manipulate matrix

            //- Check and correct zero diagonal in hole cells
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


        // Functions for imposing fringe conservation constraint

            //- Correct interpolated accepto value for fringe conservation
            void correctFringeConservation
            (
                scalarField& psiInterpolated,
                const lduMatrix& matrix
            ) const;


public:

    //- Runtime type information
    TypeName(oversetFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        oversetFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        oversetFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        oversetFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given oversetFvPatchField onto a new patch
        oversetFvPatchField
        (
            const oversetFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        oversetFvPatchField(const oversetFvPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new oversetFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        oversetFvPatchField
        (
            const oversetFvPatchField<Type>&,
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
                new oversetFvPatchField<Type>(*this, iF)
            );
        }


    // Destructor

        ~oversetFvPatchField();


    // Static member functions

        //- Find the overset patch for given field and set the coupled fringe to
        //  true (default) or false. Needed for correct fringe update of 1/A
        //  field before interpolation in the pressure equation.
        static void setCoupledFringe
        (
            GeometricField<Type, fvPatchField, volMesh>& psi,
            const bool coupled = true
        );

        //- Force the interpolation of the field, irrespective whether the
        //  coupledFringe and conservativeCorrection is specified or not
        static void oversetInterpolate
        (
            GeometricField<Type, fvPatchField, volMesh>& psi
        );


    // Member functions

        // Access

            //- Return true if coupled
            virtual bool coupled() const
            {
                return coupledFringe_;
            }

            //- Return true if the conservative correction is specified
            bool conservativeCorrection() const
            {
                return conservativeCorrection_;
            }

            //- Return neighbour field.  Dummy
            virtual tmp<Field<Type> > patchNeighbourField() const;

            //- Set coupledFringe to true or false
            void setCoupledFringe(const bool coupled)
            {
                coupledFringe_ = coupled;
            }

            //- Set conservativeCorrection to true or false
            void setConservativeCorrection(const bool conservative)
            {
                conservativeCorrection_ = conservative;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

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


        // Coupled interface functionality

            //- Transform given patch component field
            virtual void transformCoupleField
            (
                scalarField& f,
                const direction cmpt
            ) const;

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
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;

            //- Calculate patch flux
            virtual void patchFlux
            (
                GeometricField<Type, fvsPatchField, surfaceMesh>& pFlux,
                const fvMatrix<Type>& matrix
            ) const;


        //- Overset coupled interface functions


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "oversetFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
