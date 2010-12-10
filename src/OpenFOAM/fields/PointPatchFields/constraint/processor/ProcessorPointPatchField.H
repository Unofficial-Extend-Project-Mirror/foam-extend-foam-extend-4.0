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
    ProcessorPointPatchField

Description

SourceFiles
    ProcessorPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef ProcessorPointPatchField_H
#define ProcessorPointPatchField_H

#include "CoupledPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class ProcessorPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class ProcessorPointPatch,
    template<class> class MatrixType,
    class Type
>
class ProcessorPointPatchField
:
    public CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename ProcessorPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >
{
    // Private data

        //- Local reference to processor patch
        const ProcessorPointPatch& procPatch_;


        // Non-blocking parallel communications
        // Temporary: move to patch.  HJ, 15/Apr/2008

            //- Send buffer.
            //  Only sized and used when compressed or non-blocking comms used.
            mutable List<char> sendBuf_;

            //- Receive buffer.
            //  Only sized and used when compressed or non-blocking comms used.
            mutable List<char> receiveBuf_;

            //- Resize the buffer if required
            void resizeBuf(List<char>& buf, const label size) const;


    // Private member functions

        //- Raw field send function
        template<class Type2>
        void sendField
        (
            const tmp<Field<Type2> >&,
            const Pstream::commsTypes commsType = Pstream::blocking
        ) const;

        //- Raw field receive function
        template<class Type2>
        tmp<Field<Type2> > receivePointField
        (
            const Pstream::commsTypes commsType = Pstream::blocking
        ) const;

        template<class Type2>
        tmp<Field<Type2> > receiveEdgeField
        (
            const Pstream::commsTypes commsType = Pstream::blocking
        ) const;

        //- Templated add diagonal/source functions
        template<class Type2>
        void initAddFieldTempl
        (
            const Pstream::commsTypes commsType,
            const Field<Type2>&
        ) const;

        template<class Type2>
        void addFieldTempl
        (
            const Pstream::commsTypes commsType,
            Field<Type2>&
        ) const;


public:

    //- Runtime type information
    TypeName(ProcessorPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        ProcessorPointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&
        );

        //- Construct from patch, internal field and dictionary
        ProcessorPointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        ProcessorPointPatchField
        (
            const ProcessorPointPatchField
            <
                PatchField,
                Mesh,
                PointPatch,
                ProcessorPointPatch,
                MatrixType,
                Type
            >&,
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const PointPatchFieldMapper&
        );

        //- Construct as copy
        ProcessorPointPatchField
        (
            const ProcessorPointPatchField
            <
                PatchField,
                Mesh,
                PointPatch,
                ProcessorPointPatch,
                MatrixType,
                Type
            >&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const
        {
            return autoPtr<PatchField<Type> >
            (
                new ProcessorPointPatchField
                <
                    PatchField,
                    Mesh,
                    PointPatch,
                    ProcessorPointPatch,
                    MatrixType,
                    Type
                >(*this)
                );
        }

        //- Construct as copy setting internal field reference
        ProcessorPointPatchField
        (
            const ProcessorPointPatchField
            <
                PatchField,
                Mesh,
                PointPatch,
                ProcessorPointPatch,
                MatrixType,
                Type
            >&,
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
                new ProcessorPointPatchField
                <
                    PatchField,
                    Mesh,
                    PointPatch,
                    ProcessorPointPatch,
                    MatrixType,
                    Type
                >(*this, iF)
            );
        }


    // Destructor

        virtual ~ProcessorPointPatchField();


    // Member functions

        // Access

            //- Return true if running parallel
            virtual bool coupled() const
            {
                if (Pstream::parRun())
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }

            //- Does the patch field perform the transfromation
            virtual bool doTransform() const
            {
                return
                   !(
                        procPatch_.procPolyPatch().parallel()
                     || pTraits<Type>::rank == 0
                    );
            }

            //- Initialise evaluation of the patch field
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );

            //- Init add field to internal field
            virtual void initAddField() const;

            //- Add field to internal field
            virtual void addField(Field<Type>&) const;

            //- Set boundary condition to matrix.  Used to save the processor
            //  patch off-diagonal coefficients
            virtual void setBoundaryCondition
            (
                Map<typename MatrixType<Type>::ConstraintType> &
            ) const;


        // Matrix construction (completion) functionality

            //- Add diagonal/source contributions
            virtual void initAddDiag(const scalarField&) const;
            virtual void initAddSource(const scalarField&) const;

            virtual void addDiag(scalarField&) const;
            virtual void addSource(scalarField&) const;

            //- Add upper/lower contributions
            virtual void initAddUpperLower(const scalarField&) const;

            virtual void addUpperLower(scalarField&) const;

            virtual void eliminateUpperLower(scalarField&) const;

            //- Get the cut edge coefficients in Amul order
            virtual tmp<scalarField> cutBouCoeffs(const lduMatrix&) const;

            //- Get the cut edge coefficients in Tmul order
            virtual tmp<scalarField> cutIntCoeffs(const lduMatrix&) const;


        // Coupled interface functionality

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "ProcessorPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
