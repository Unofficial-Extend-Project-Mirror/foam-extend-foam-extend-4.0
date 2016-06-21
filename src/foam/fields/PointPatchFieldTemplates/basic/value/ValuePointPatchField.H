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
    ValuePointPatchField

Description

SourceFiles
    ValuePointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef ValuePointPatchField_H
#define ValuePointPatchField_H

#include "PointPatchFieldTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class ValuePointPatchField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
class ValuePointPatchField
:
    public PatchField<Type>,
    public Field<Type>
{
    // Private member functions

        void checkFieldSize() const;


public:

    //- Runtime type information
    TypeName("value");


    // Constructors

        //- Construct from patch and internal field
        ValuePointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&
        );

        //- Construct from patch, internal field and dictionary
        ValuePointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        ValuePointPatchField
        (
            const ValuePointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>&,
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const PointPatchFieldMapper&
        );

        //- Construct as copy
        ValuePointPatchField
        (
            const ValuePointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const
        {
            return autoPtr<PatchField<Type> >
            (
                new ValuePointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        ValuePointPatchField
        (
            const ValuePointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>&,
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
                new ValuePointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return size
            label size() const
            {
                return Field<Type>::size();
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const PointPatchFieldMapper&
            );

            //- Reverse map the given PointPatchField onto
            //  this PointPatchField
            virtual void rmap
            (
                const PointPatchField
                    <PatchField, Mesh, PointPatch, MatrixType, Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Does it store field data
            virtual bool storesFieldData() const
            {
                return true;
            }

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );

            //- Update boundary value
            virtual void updateBoundaryField()
            {}

            //- Set boundary condition to matrix
            virtual void setBoundaryCondition
            (
                Map<typename MatrixType<Type>::ConstraintType> &
            ) const
            {}


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        //- Assignment operators

            virtual void operator=
            (
                const ValuePointPatchField
                    <PatchField, Mesh, PointPatch, MatrixType, Type>&
            );

            virtual void operator=
            (
                const PointPatchField
                    <PatchField, Mesh, PointPatch, MatrixType, Type>&
            );

            virtual void operator=(const Field<Type>&);
            virtual void operator=(const Type&);


        // Force an assignment irrespective of form of patch

            virtual void operator==
            (
                const ValuePointPatchField
                    <PatchField, Mesh, PointPatch, MatrixType, Type>&
            );

            virtual void operator==
            (
                const PointPatchField
                    <PatchField, Mesh, PointPatch, MatrixType, Type>&
            );

            virtual void operator==(const Field<Type>&);
            virtual void operator==(const Type&);
};


// This function is added to overide the hack in PointPatchField.H which
//  enables simple backward compatability with versions using referenceLevel
// in GeometricField
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
tmp<Field<Type> > operator+
(
    const ValuePointPatchField
    <PatchField, Mesh, PointPatch, MatrixType, Type>& vsppf,
    const Type& t
)
{
    return static_cast<const Field<Type>&>(vsppf) + t;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "ValuePointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
