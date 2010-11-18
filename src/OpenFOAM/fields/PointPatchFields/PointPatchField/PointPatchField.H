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
    PointPatchField

Description
    Abstract base class for point-mesh patch fields.

    The base-field does not store values as they are part of the
    "internal field".  There are derived classes to store constraint values
    e.g. FixedValuePointPatchField derived from the generic
    ValuePointPatchField which ensures the values in the "internal field"
    are reset to the fixed-values by applying the stored values.

SourceFiles
    PointPatchField.C
    newPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef PointPatchField_H
#define PointPatchField_H

#include "DimensionedField.H"
#include "fieldTypes.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
class PointPatchField;

template<class T> class Map;

class objectRegistry;
class dictionary;
class PointPatchFieldMapper;
class lduMatrix;

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
    >
Ostream& operator<<
(
    Ostream&,
    const PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>&
);


/*---------------------------------------------------------------------------*\
                        Class PointPatchField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
class PointPatchField
{

public:

    typedef PointPatch Patch;
    typedef PatchField<Type> pointPatchTypeField;


private:

    // Private data

        //- Reference to patch
        const PointPatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, Mesh>& internalField_;

        //- Update index used so that updateCoeffs is called only once
        bool updated_;


public:

    //- Runtime type information
    TypeName("PointPatchField");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointPatchTypeField,
            PointPatch,
            (
                const PointPatch& p,
                const DimensionedField<Type, Mesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointPatchTypeField,
            patchMapper,
            (
                const PointPatchField
                    <PatchField, Mesh, PointPatch, MatrixType, Type>& ptf,
                const PointPatch& p,
                const DimensionedField<Type, Mesh>& iF,
                const PointPatchFieldMapper& m
            ),
            (dynamic_cast<const pointPatchTypeFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointPatchTypeField,
            dictionary,
            (
                const PointPatch& p,
                const DimensionedField<Type, Mesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        PointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&
        );

        //- Construct as copy
        PointPatchField
        (
            const PointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const = 0;

        //- Construct as copy setting internal field reference
        PointPatchField
        (
            const PointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>&,
            const DimensionedField<Type, Mesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<PatchField<Type> > clone
        (
            const DimensionedField<Type, Mesh>& iF
        ) const = 0;


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static autoPtr<PatchField<Type> > New
        (
            const word&,
            const PointPatch&,
            const DimensionedField<Type, Mesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given PointPatchField mapped onto a new patch
        static autoPtr<PatchField<Type> > New
        (
            const PointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>&,
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const PointPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static autoPtr<PatchField<Type> > New
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const dictionary&
        );


        //- Return a pointer to a new CalculatedPointPatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static
        autoPtr
        <
            PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>
        >
        NewCalculatedType
        (
            const PointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type2>&
        );


    // Destructor

        virtual ~PointPatchField()
        {}


    // Member functions

        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return size
            label size() const
            {
                return patch().size();
            }

            //- Return patch
            const PointPatch& patch() const
            {
                return patch_;
            }

            //- Return internal field reference
            const DimensionedField<Type, Mesh>& dimensionedInternalField() const
            {
                return internalField_;
            }

            //- Return internal field reference
            const Field<Type>& internalField() const
            {
                return internalField_;
            }

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return false;
            }

            //- Return true if the boundary condition has already been updated
            bool updated() const
            {
                return updated_;
            }

            //- Return field created from appropriate internal field values
            tmp<Field<Type> > patchInternalField() const;

            //- Return field created from appropriate internal field values
            //  given internal field reference
            template<class Type1>
            tmp<Field<Type1> > patchInternalField
            (
                const Field<Type1>& iF
            ) const;

            //- Given the internal field and a patch field,
            //  add the patch field to the internal field
            template<class Type1>
            void addToInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF
            ) const;

            //- Given the internal field and a patch field,
            //  set the patch field in the internal field
            template<class Type1>
            void setInInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF
            ) const;

            //- Return the type of the calculated for of PointPatchField
            static const word& calculatedType();


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const PointPatchFieldMapper&
            )
            {}

            //- Reverse map the given PointPatchField onto
            //  this PointPatchField
            virtual void rmap
            (
                const PointPatchField
                    <PatchField, Mesh, PointPatch, MatrixType, Type>&,
                const labelList&
            )
            {}


        // Evaluation functions

            //- Does it store field data
            virtual bool storesFieldData() const
            {
                return false;
            }

            //- Does this patchField correspond to a pointTypeField
            bool isPointField() const;

            //- Check that this patchField corresponds to a pointTypeField,
            // if not abort!
            void checkPointField() const;


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs()
            {
                updated_ = true;
            }

            //- Initialise evaluation of the patch field (do nothing)
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            )
            {}

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );

            //- Update boundary value
            virtual void updateBoundaryField()
            {}

            //- Init add field to internal field
            virtual void initAddField() const
            {}

            //- Add field to internal field
            virtual void addField(Field<Type>&) const
            {
                notImplemented(type() + "::addField(Field<Type>&) const");
            }

            //- Add boundary source for gradient-type conditions
            virtual void addBoundarySourceDiag
            (
                MatrixType<Type>&
            ) const
            {}

            //- Set boundary condition to matrix
            virtual void setBoundaryCondition
            (
                Map<typename MatrixType<Type>::ConstraintType>&
            ) const
            {}


        // Matrix construction (completion) functionality

            //- Add diagonal/source contributions
            virtual void initAddDiag(const scalarField&) const
            {}

            virtual void initAddSource(const scalarField&) const
            {}

            virtual void addDiag(scalarField&) const
            {
                notImplemented(type() + "::addDiag(scalarField&) const");
            }

            virtual void addSource(scalarField&) const
            {
                notImplemented(type() + "::addSource(scalarField&) const");
            }

            //- Get the cut edge coefficients in Amul order
            virtual tmp<scalarField> cutBouCoeffs(const lduMatrix&) const
            {
                return tmp<scalarField>(new scalarField(this->size(), 0));
            }

            //- Get the cut edge coefficients in Tmul order
            virtual tmp<scalarField> cutIntCoeffs(const lduMatrix&) const
            {
                return tmp<scalarField>(new scalarField(this->size(), 0));
            }

            //- Add upper/lower contributions
            virtual void initAddUpperLower(const scalarField&) const
            {}

            virtual void addUpperLower(scalarField&) const
            {
                notImplemented(type() + "::addUpperLower(scalarField&) const");
            }

            virtual void eliminateUpperLower(scalarField&) const
            {
                notImplemented
                (
                    type()
                  + "::eliminateUpperLower(scalarField&) const"
                );
            }


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=
        (
            const PointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>&
        )
        {}

        virtual void operator+=
        (
            const PointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>&
        )
        {}

        virtual void operator-=
        (
            const PointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>&
        )
        {}

        virtual void operator*=
        (
            const PointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, scalar>&
        )
        {}

        virtual void operator/=
        (
            const PointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, scalar>&
        )
        {}

        virtual void operator=(const Field<Type>&){}
        virtual void operator+=(const Field<Type>&){}
        virtual void operator-=(const Field<Type>&){}

        virtual void operator*=(const Field<scalar>&){}
        virtual void operator/=(const Field<scalar>&){}

        virtual void operator=(const Type&){}
        virtual void operator+=(const Type&){}
        virtual void operator-=(const Type&){}
        virtual void operator*=(const scalar){}
        virtual void operator/=(const scalar){}


        // Force an assignment irrespective of form of patch
        // By default these do nothing unless the patch actually has boundary
        // values

        virtual void operator==
        (
            const PointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>&
        )
        {}

        virtual void operator==(const Field<Type>&){}
        virtual void operator==(const Type&){}


    // Ostream operator

        friend Ostream& operator<<
        <PatchField, Mesh, PointPatch, MatrixType, Type>
        (
            Ostream&,
            const PointPatchField
                <PatchField, Mesh, PointPatch, MatrixType, Type>&
        );
};


// This function is added as a hack to enable simple backward compatability
// with verions using referenceLevel in GeometicField
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    template<class> class MatrixType,
    class Type
>
const PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>& operator+
(
    const PointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>& ppf,
    const Type&
)
{
    return ppf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PointPatchFieldFunctions.H"

#ifdef NoRepository
#   include "PointPatchField.C"
#   include "CalculatedPointPatchField.H"
#endif


#define makePointPatchTypeFieldTypeName(typePatchTypeField)                   \
                                                                              \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);

#define makePointPatchFieldsTypeName(typePatchField, pointPatchType)          \
                                                                              \
makePointPatchTypeFieldTypeName(typePatchField##pointPatchType##ScalarField); \
makePointPatchTypeFieldTypeName(typePatchField##pointPatchType##VectorField); \
makePointPatchTypeFieldTypeName(typePatchField##pointPatchType##SphericalTensorField); \
makePointPatchTypeFieldTypeName(typePatchField##pointPatchType##SymmTensorField); \
makePointPatchTypeFieldTypeName(typePatchField##pointPatchType##TensorField);

#define makePointPatchTypeField(PatchTypeField, typePatchTypeField)           \
                                                                              \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);                   \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    PatchTypeField, typePatchTypeField, PointPatch                            \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    PatchTypeField,                                                           \
    typePatchTypeField,                                                       \
    patchMapper                                                               \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    PatchTypeField, typePatchTypeField, dictionary                            \
);


#define makePointPatchFields(type, pointPatch, pointPatchType)                \
                                                                              \
makePointPatchTypeField                                                       \
(                                                                             \
    pointPatch##ScalarField,                                                  \
    type##pointPatchType##ScalarField                                         \
);                                                                            \
                                                                              \
makePointPatchTypeField                                                       \
(                                                                             \
    pointPatch##VectorField,                                                  \
    type##pointPatchType##VectorField                                         \
);                                                                            \
                                                                              \
makePointPatchTypeField                                                       \
(                                                                             \
    pointPatch##SphericalTensorField,                                         \
    type##pointPatchType##SphericalTensorField                                \
);                                                                            \
                                                                              \
makePointPatchTypeField                                                       \
(                                                                             \
    pointPatch##SymmTensorField,                                              \
    type##pointPatchType##SymmTensorField                                     \
);                                                                            \
                                                                              \
makePointPatchTypeField                                                       \
(                                                                             \
    pointPatch##TensorField,                                                  \
    type##pointPatchType##TensorField                                         \
);
             

#define makePointPatchFieldTypedefs(type, Type, mesh, pointPatch, patchPatch, matrixType, pointPatchType) \
                                                                              \
typedef Type##PointPatchField                                                 \
    <pointPatch##Field, mesh, patchPatch, matrixType, scalar>                 \
    type##pointPatchType##ScalarField;                                        \
typedef Type##PointPatchField                                                 \
    <pointPatch##Field, mesh, patchPatch, matrixType, vector>                 \
    type##pointPatchType##VectorField;                                        \
typedef Type##PointPatchField                                                 \
    <pointPatch##Field, mesh, patchPatch, matrixType, sphericalTensor>        \
    type##pointPatchType##SphericalTensorField;                               \
typedef Type##PointPatchField                                                 \
    <pointPatch##Field, mesh, patchPatch, matrixType, symmTensor>             \
    type##pointPatchType##SymmTensorField;                                    \
typedef Type##PointPatchField                                                 \
    <pointPatch##Field, mesh, patchPatch, matrixType, tensor>                 \
    type##pointPatchType##TensorField;

#define makeTypePointPatchFieldTypedefs(type, Type, mesh, pointPatch, patchPatch, pointPatchType, matrixType, typePointPatch) \
                                                                              \
typedef Type##PointPatchField                                                 \
    <pointPatch##Field, mesh, patchPatch, typePointPatch, matrixType, scalar> \
    type##pointPatchType##ScalarField;                                        \
typedef Type##PointPatchField                                                 \
    <pointPatch##Field, mesh, patchPatch, typePointPatch, matrixType, vector> \
    type##pointPatchType##VectorField;                                        \
typedef Type##PointPatchField                                                 \
    <pointPatch##Field, mesh, patchPatch, typePointPatch, matrixType, sphericalTensor>\
    type##pointPatchType##SphericalTensorField;                               \
typedef Type##PointPatchField                                                 \
    <pointPatch##Field, mesh, patchPatch, typePointPatch, matrixType, symmTensor>\
    type##pointPatchType##SymmTensorField;                                    \
typedef Type##PointPatchField                                                 \
    <pointPatch##Field, mesh, patchPatch, typePointPatch, matrixType, tensor> \
    type##pointPatchType##TensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
