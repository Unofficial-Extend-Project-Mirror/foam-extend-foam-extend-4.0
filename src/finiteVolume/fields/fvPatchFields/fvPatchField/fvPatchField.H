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
    Foam::fvPatchField

Description
    Abstract base class with a fat-interface to all derived classes
    covering all possible ways in which they might be used.

    The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.

    The next level of derivation covers all the specialised types with
    specific evaluation proceedures, particularly with respect to specific
    fields.

SourceFiles
    fvPatchField.C
    newFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fvPatchField_H
#define fvPatchField_H

#include "fvPatch.H"
#include "DimensionedField.H"
#include "debugSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class objectRegistry;
class dictionary;
class fvPatchFieldMapper;
class volMesh;
class surfaceMesh;

// Forward declaration of friend functions and operators

template<class Type>
class fvPatchField;

template<class Type>
class fvsPatchField;

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

template<class Type>
class fvMatrix;

template<class Type>
Ostream& operator<<(Ostream&, const fvPatchField<Type>&);


/*---------------------------------------------------------------------------*\
                        Class fvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvPatchField
:
    public Field<Type>
{
    // Private data

        //- Reference to patch
        const fvPatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, volMesh>& internalField_;

        //- Update index used so that updateCoeffs is called only once during
        //  the construction of the matrix
        bool updated_;

        //- Optional patch type, used to allow specified boundary conditions
        //  to be applied to constraint patches by providing the constraint
        //  patch type as 'patchType'
        word patchType_;


public:

    typedef fvPatch Patch;

    //- Runtime type information
    TypeName("fvPatchField");

    //- Debug switch to disallow the use of genericFvPatchField
    static debug::debugSwitch disallowGenericFvPatchField;


    // Declare run-time constructor selection tables

#ifndef SWIG
        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            patch,
            (
                const fvPatch& p,
                const DimensionedField<Type, volMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            patchMapper,
            (
                const fvPatchField<Type>& ptf,
                const fvPatch& p,
                const DimensionedField<Type, volMesh>& iF,
                const fvPatchFieldMapper& m
            ),
            (dynamic_cast<const fvPatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            dictionary,
            (
                const fvPatch& p,
                const DimensionedField<Type, volMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );
#endif

    // Constructors

        //- Construct from patch and internal field
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&,
            const bool valueRequired = false
        );

        //- Construct by mapping the given fvPatchField onto a new patch
        fvPatchField
        (
            const fvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fvPatchField(const fvPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >(new fvPatchField<Type>(*this));
        }

        //- Construct as copy setting internal field reference
        fvPatchField
        (
            const fvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >(new fvPatchField<Type>(*this, iF));
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<fvPatchField<Type> > New
        (
            const word&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given fvPatchField mapped onto a new patch
        static tmp<fvPatchField<Type> > New
        (
            const fvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<fvPatchField<Type> > New
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedFvPatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static tmp<fvPatchField<Type> > NewCalculatedType
        (
            const fvPatchField<Type2>&
        );


    // Destructor

        virtual ~fvPatchField<Type>()
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
            const DimensionedField<Type, volMesh>&
            dimensionedInternalField() const
            {
                return internalField_;
            }

            //- Return internal field reference
            const Field<Type>& internalField() const
            {
                return internalField_;
            }

            //- Return the type of the calculated for type of fvPatchField
            static const word& calculatedType();

            //- Return patch type (optional on read)
            const word& patchType() const
            {
                return patchType_;
            }

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

            //- Return true if the boundary condition has already been updated
            bool updated() const
            {
                return updated_;
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

            //- Return patch-normal gradient
            virtual tmp<Field<Type> > snGrad() const;

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs()
            {
                updated_ = true;
            }

            //- Return internal field next to patch as patch field
            virtual tmp<Field<Type> > patchInternalField() const;

            //- Return patchField on the opposite patch of a coupled patch
            virtual tmp<Field<Type> > patchNeighbourField() const
            {
                notImplemented(type() + "patchNeighbourField()");
                return *this;
            }

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
            );


            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueInternalCoeffs
            (
                const tmp<Field<scalar> >&
            ) const
            {
                notImplemented
                (
                    type()
                  + "::valueInternalCoeffs(const tmp<Field<scalar> >&)"
                );
                return *this;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueBoundaryCoeffs
            (
                const tmp<Field<scalar> >&
            ) const
            {
                notImplemented
                (
                    type()
                  + "::valueBoundaryCoeffs(const tmp<Field<scalar> >&)"
                );
                return *this;
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientInternalCoeffs() const
            {
                notImplemented(type() + "::gradientInternalCoeffs()");
                return *this;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientBoundaryCoeffs() const
            {
                notImplemented(type() + "::gradientBoundaryCoeffs()");
                return *this;
            }


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

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<Type>& matrix);


        // I-O

            //- Write
            virtual void write(Ostream&) const;

            //- Helper function to write the keyword and entry only if the
            //  values are not equal. The value is then output as value2
            template<class EntryType>
            void writeEntryIfDifferent
            (
                Ostream& os,
                const word& entryName,
                const EntryType& value1,
                const EntryType& value2
            ) const;


        // Check

            //- Check fvPatchField<Type> against given fvPatchField<Type>
            void check(const fvPatchField<Type>&) const;


    // Member operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const fvPatchField<Type>&);
        virtual void operator+=(const fvPatchField<Type>&);
        virtual void operator-=(const fvPatchField<Type>&);
        virtual void operator*=(const fvPatchField<scalar>&);
        virtual void operator/=(const fvPatchField<scalar>&);

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

        virtual void operator==(const fvPatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);


    // Ostream operator

#ifndef SWIG
        friend Ostream& operator<< <Type>(Ostream&, const fvPatchField<Type>&);
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvPatchField.C"
#   include "calculatedFvPatchField.H"
#endif


#define addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)   \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    PatchTypeField, typePatchTypeField, patch                                 \
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


#define makePatchTypeFieldTypeName(type)                                      \
                                                                              \
defineNamedTemplateTypeNameAndDebug(type, 0);


#define makePatchFieldsTypeName(type)                                         \
                                                                              \
makePatchTypeFieldTypeName(type##FvPatchScalarField);                         \
makePatchTypeFieldTypeName(type##FvPatchVectorField);                         \
makePatchTypeFieldTypeName(type##FvPatchSphericalTensorField);                \
makePatchTypeFieldTypeName(type##FvPatchSymmTensorField);                     \
makePatchTypeFieldTypeName(type##FvPatchSymmTensor4thOrderField);             \
makePatchTypeFieldTypeName(type##FvPatchDiagTensorField);                     \
makePatchTypeFieldTypeName(type##FvPatchTensorField);


#define makePatchTypeField(PatchTypeField, typePatchTypeField)                \
                                                                              \
defineTypeNameAndDebug(typePatchTypeField, 0);                                \
                                                                              \
addToPatchFieldRunTimeSelection                                               \
(                                                                             \
    PatchTypeField, typePatchTypeField                                        \
);


#define makeTemplatePatchTypeField(PatchTypeField, typePatchTypeField)        \
                                                                              \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);                   \
                                                                              \
addToPatchFieldRunTimeSelection                                               \
(                                                                             \
    PatchTypeField, typePatchTypeField                                        \
);


#define makePatchFields(type)                                                 \
                                                                              \
makeTemplatePatchTypeField(fvPatchScalarField, type##FvPatchScalarField);     \
makeTemplatePatchTypeField(fvPatchVectorField, type##FvPatchVectorField);     \
makeTemplatePatchTypeField                                                    \
(                                                                             \
    fvPatchSphericalTensorField,                                              \
    type##FvPatchSphericalTensorField                                         \
);                                                                            \
makeTemplatePatchTypeField                                                    \
(                                                                             \
    fvPatchSymmTensorField,                                                   \
    type##FvPatchSymmTensorField                                              \
);                                                                            \
makeTemplatePatchTypeField                                                    \
(                                                                             \
    fvPatchSymmTensor4thOrderField,                                           \
    type##FvPatchSymmTensor4thOrderField                                      \
);                                                                            \
makeTemplatePatchTypeField                                                    \
 (                                                                            \
    fvPatchDiagTensorField,                                                   \
    type##FvPatchDiagTensorField                                              \
);                                                                            \
makeTemplatePatchTypeField(fvPatchTensorField, type##FvPatchTensorField);


#define makePatchTypeFieldTypedefs(type)                                      \
                                                                              \
typedef type##FvPatchField<scalar> type##FvPatchScalarField;                  \
typedef type##FvPatchField<vector> type##FvPatchVectorField;                  \
typedef type##FvPatchField<sphericalTensor>                                   \
    type##FvPatchSphericalTensorField;                                        \
typedef type##FvPatchField<symmTensor> type##FvPatchSymmTensorField;          \
typedef type##FvPatchField<symmTensor4thOrder>                                \
    type##FvPatchSymmTensor4thOrderField;                                     \
typedef type##FvPatchField<diagTensor>                                        \
    type##FvPatchDiagTensorField;                                             \
typedef type##FvPatchField<tensor> type##FvPatchTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
