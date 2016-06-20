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
    GlobalPointPatchField

Description
    Patch field containing globally shared points in parallel decomposition

SourceFiles
    GlobalPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef GlobalPointPatchField_H
#define GlobalPointPatchField_H

#include "CoupledPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class GlobalPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class GlobalPointPatch,
    template<class> class MatrixType,
    class Type
>
class GlobalPointPatchField
:
    public CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename GlobalPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >
{
    // Private data

        //- Local reference to processorPoint patch
        const GlobalPointPatch& globalPointPatch_;


    // Private member functions

        //- Reduce the field and extract the local values
        template<class Type2>
        tmp<Field<Type2> > reduceExtractPoint(const tmp<Field<Type2> >&) const;

        template<class Type2>
        tmp<Field<Type2> > reduceExtractEdge(const tmp<Field<Type2> >&) const;

        template<class Type2>
        void addFieldTempl(Field<Type2>&) const;


public:

    //- Runtime type information
    TypeName(GlobalPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        GlobalPointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&
        );

        //- Construct from patch, internal field and dictionary
        GlobalPointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        GlobalPointPatchField
        (
            const GlobalPointPatchField
            <PatchField, Mesh, PointPatch, GlobalPointPatch, MatrixType, Type>&,
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const PointPatchFieldMapper&
        );

        //- Construct as copy
        GlobalPointPatchField
        (
            const GlobalPointPatchField
            <PatchField, Mesh, PointPatch, GlobalPointPatch, MatrixType, Type>&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const
        {
            return autoPtr<PatchField<Type> >
            (
                new GlobalPointPatchField
                <
                    PatchField,
                    Mesh,
                    PointPatch,
                    GlobalPointPatch,
                    MatrixType,
                    Type
                >(*this)
            );
        }

        //- Construct as copy setting internal field reference
        GlobalPointPatchField
        (
            const GlobalPointPatchField
            <PatchField, Mesh, PointPatch, GlobalPointPatch, MatrixType, Type>&,
            const DimensionedField<Type, Mesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<PatchField<Type> > clone
        (
            const DimensionedField<Type, Mesh>& iF
        ) const
        {
            return autoPtr<PatchField<Type> >
            (
                new GlobalPointPatchField
                <
                    PatchField,
                    Mesh,
                    PointPatch,
                    GlobalPointPatch,
                    MatrixType,
                    Type
                >(*this, iF)
            );
        }


    // Destructor

        ~GlobalPointPatchField();


    // Member functions

        // Access

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );

            //- Add field to internal field
            virtual void addField(Field<Type>&) const;

            //- Set boundary condition to matrix.  Used to save the processor
            //  patch off-diagonal coefficients
            virtual void setBoundaryCondition
            (
                Map<typename MatrixType<Type>::ConstraintType> &
            ) const;


        // Matrix construction (completion) functionality

            virtual void addDiag(scalarField&) const;
            virtual void addSource(scalarField&) const;

            //- Add upper/lower contributions
            virtual void addUpperLower(scalarField&) const;

            virtual void eliminateUpperLower(scalarField&) const;

            //- Get the cut edge coefficients in Amul order
            virtual tmp<scalarField> cutBouCoeffs(const lduMatrix&) const;

            //- Get the cut edge coefficients in Tmul order
            virtual tmp<scalarField> cutIntCoeffs(const lduMatrix&) const;


        // Coupled interface functionality

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

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "GlobalPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
