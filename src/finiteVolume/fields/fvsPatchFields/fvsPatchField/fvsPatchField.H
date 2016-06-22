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
    Foam::fvsPatchField

Description
    An abstract base class with a fat-interface to all derived classes
    covering all possible ways in which they might be used.

    The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.

    The next level of derivation covers all the specialised typed with
    specific evaluation proceedures, particularly with respect to specific
    fields.

SourceFiles
    fvsPatchField.C
    newFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fvsPatchField_H
#define fvsPatchField_H

#include "fvPatch.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class objectRegistry;
class dictionary;
class fvPatchFieldMapper;
class surfaceMesh;

template<class Type>
class fvPatchField;

// Forward declaration of friend functions and operators

template<class Type>
class fvsPatchField;

template<class Type>
Ostream& operator<<(Ostream&, const fvsPatchField<Type>&);


/*---------------------------------------------------------------------------*\
                           Class patch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvsPatchField
:
    public Field<Type>
{
    // Private data

        //- Reference to patch
        const fvPatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, surfaceMesh>& internalField_;


public:

    typedef fvPatch Patch;


    //- Runtime type information
    TypeName("fvsPatchField");

    //- Debug switch to disallow the use of
    static debug::debugSwitch disallowDefaultFvsPatchField;


    // Declare run-time constructor selection tables

#ifndef SWIG
        declareRunTimeSelectionTable
        (
            tmp,
            fvsPatchField,
            patch,
            (
                const fvPatch& p,
                const DimensionedField<Type, surfaceMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            fvsPatchField,
            patchMapper,
            (
                const fvsPatchField<Type>& ptf,
                const fvPatch& p,
                const DimensionedField<Type, surfaceMesh>& iF,
                const fvPatchFieldMapper& m
            ),
            (dynamic_cast<const fvsPatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            fvsPatchField,
            dictionary,
            (
                const fvPatch& p,
                const DimensionedField<Type, surfaceMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );
#endif

    // Constructors

        //- Construct from patch and internal field
        fvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch and internal field and patch field
        fvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        fvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&,
            const bool valueRequired = false
        );

        //- Construct by mapping the given fvsPatchField onto a new patch
        fvsPatchField
        (
            const fvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fvsPatchField(const fvsPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type> > clone() const
        {
            return tmp<fvsPatchField<Type> >(new fvsPatchField<Type>(*this));
        }

        //- Construct as copy setting internal field reference
        fvsPatchField
        (
            const fvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type> > clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type> >
            (
                new fvsPatchField<Type>(*this, iF)
            );
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<fvsPatchField<Type> > New
        (
            const word&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given fvsPatchField mapped onto a new patch
        static tmp<fvsPatchField<Type> > New
        (
            const fvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<fvsPatchField<Type> > New
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedFvsPatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static tmp<fvsPatchField<Type> > NewCalculatedType
        (
            const fvsPatchField<Type2>&
        );


    // Destructor

        virtual ~fvsPatchField<Type>()
        {}


    // Member functions

        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Lookup and return the patchField of the named field from the
            //  field's objectRegistry.
            //  N.B.  The dummy pointer arguments are used if this function is
            //  instantiated within a templated function to avoid a bug in gcc.
            //  See inletOutletFvPatchField.C and outletInletFvPatchField.C
            template<class GeometricField, class Type2>
            const typename GeometricField::PatchFieldType& lookupPatchField
            (
                const word& name,
                const GeometricField* = NULL,
                const Type2* = NULL
            ) const;

            //- Return patch
            const fvPatch& patch() const
            {
                return patch_;
            }

            //- Return dimensioned internal field reference
            const DimensionedField<Type, surfaceMesh>&
            dimensionedInternalField() const
            {
                return internalField_;
            }

            //- Return internal field reference
            const Field<Type>& internalField() const
            {
                return internalField_;
            }

            //- Return the type of the calculated for of fvsPatchField
            static const word& calculatedType();

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return false;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvsPatchField onto this fvsPatchField
            virtual void rmap
            (
                const fvsPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            )
            {}

            //- Evaluate the patch field, sets Updated to false
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            )
            {}


        //- Write
        virtual void write(Ostream&) const;


        // Check

            //- Check fvsPatchField<Type> against given fvsPatchField<Type>
            void check(const fvsPatchField<Type>&) const;

            //- Check fvsPatchField<Type> against given fvPatchField<Type>
            void check(const fvPatchField<Type>&) const;


    // Member operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const fvsPatchField<Type>&);
        virtual void operator+=(const fvsPatchField<Type>&);
        virtual void operator-=(const fvsPatchField<Type>&);
        virtual void operator*=(const fvsPatchField<scalar>&);
        virtual void operator/=(const fvsPatchField<scalar>&);

        virtual void operator=(const fvPatchField<Type>&);

        virtual void operator+=(const Field<Type>&);
        virtual void operator-=(const Field<Type>&);

        virtual void operator*=(const Field<scalar>&);
        virtual void operator/=(const Field<scalar>&);

        virtual void operator=(const Type&);
        virtual void operator+=(const Type&);
        virtual void operator-=(const Type&);
        virtual void operator*=(const scalar);
        virtual void operator/=(const scalar);


        // Force an assignment irrespective of form of patch

        virtual void operator==(const fvsPatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);


    // Ostream operator

#ifndef SWIG
        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const fvsPatchField<Type>&
        );
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvsPatchField.C"
#   include "calculatedFvsPatchField.H"
#endif


#define makeFvsPatchTypeFieldTypeName(type)                                \
                                                                           \
defineNamedTemplateTypeNameAndDebug(type, 0);

#define makeFvsPatchFieldsTypeName(type)                                   \
                                                                           \
makeFvsPatchTypeFieldTypeName(type##FvsPatchScalarField);                  \
makeFvsPatchTypeFieldTypeName(type##FvsPatchVectorField);                  \
makeFvsPatchTypeFieldTypeName(type##FvsPatchSphericalTensorField);         \
makeFvsPatchTypeFieldTypeName(type##FvsPatchSymmTensorField);              \
makeFvsPatchTypeFieldTypeName(type##FvsPatchSymmTensor4thOrderField);      \
makeFvsPatchTypeFieldTypeName(type##FvsPatchDiagTensorField);              \
makeFvsPatchTypeFieldTypeName(type##FvsPatchTensorField);

#define makeFvsPatchTypeField(PatchTypeField, typePatchTypeField)          \
                                                                           \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);                \
                                                                           \
addToRunTimeSelectionTable                                                 \
(                                                                          \
    PatchTypeField, typePatchTypeField, patch                              \
);                                                                         \
                                                                           \
addToRunTimeSelectionTable                                                 \
(                                                                          \
    PatchTypeField,                                                        \
    typePatchTypeField,                                                    \
    patchMapper                                                            \
);                                                                         \
                                                                           \
addToRunTimeSelectionTable                                                 \
(                                                                          \
    PatchTypeField, typePatchTypeField, dictionary                         \
);


#define makeFvsPatchFields(type)                                           \
                                                                           \
makeFvsPatchTypeField(fvsPatchScalarField, type##FvsPatchScalarField);     \
makeFvsPatchTypeField(fvsPatchVectorField, type##FvsPatchVectorField);     \
makeFvsPatchTypeField                                                      \
(                                                                          \
    fvsPatchSphericalTensorField,                                          \
    type##FvsPatchSphericalTensorField                                     \
);                                                                         \
makeFvsPatchTypeField(fvsPatchSymmTensorField, type##FvsPatchSymmTensorField); \
makeFvsPatchTypeField                                                      \
(                                                                          \
    fvsPatchSymmTensor4thOrderField,                                       \
    type##FvsPatchSymmTensor4thOrderField                                  \
);                                                                         \
makeFvsPatchTypeField                                                      \
(                                                                          \
    fvsPatchDiagTensorField,                                               \
    type##FvsPatchDiagTensorField                                          \
);                                                                         \
makeFvsPatchTypeField(fvsPatchTensorField, type##FvsPatchTensorField);


#define makeFvsPatchTypeFieldTypedefs(type)                                \
                                                                           \
typedef type##FvsPatchField<scalar> type##FvsPatchScalarField;             \
typedef type##FvsPatchField<vector> type##FvsPatchVectorField;             \
typedef type##FvsPatchField<sphericalTensor>                               \
    type##FvsPatchSphericalTensorField;                                    \
typedef type##FvsPatchField<symmTensor> type##FvsPatchSymmTensorField;     \
typedef type##FvsPatchField<symmTensor4thOrder>                            \
    type##FvsPatchSymmTensor4thOrderField;                                 \
typedef type##FvsPatchField<diagTensor>                                    \
    type##FvsPatchDiagTensorField;                                         \
typedef type##FvsPatchField<tensor> type##FvsPatchTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
