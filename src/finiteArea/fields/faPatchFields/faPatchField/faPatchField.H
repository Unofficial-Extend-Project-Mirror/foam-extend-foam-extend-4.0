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
    faPatchField

Description
    faPatchField<Type> abstract base class.  This class gives a fat-interface
    to all derived classes covering all possible ways in which they might be
    used.  The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.  The next
    level of derivation covers all the specialised typed with specific
    evaluation proceedures, particularly with respect to specific fields.

SourceFiles
    faPatchField.C
    newPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef faPatchField_H
#define faPatchField_H

#include "faPatch.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

class objectRegistry;
class dictionary;
class faPatchFieldMapper;
class areaMesh;


// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class faPatchField;

template<class Type>
Ostream& operator<<(Ostream&, const faPatchField<Type>&);

/*---------------------------------------------------------------------------*\
                         Class faPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class faPatchField
:
    public Field<Type>
{
    // Private data

        //- Reference to a patch
        const faPatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, areaMesh>& internalField_;

        //- Update index used so that updateCoeffs is called only once during
        //  the construction of the matrix
        bool updated_;


public:

    typedef faPatch Patch;


    //- Runtime type information
    TypeName("faPatchField");

    //- Debug switch to disallow the use of 
    static int disallowDefaultFaPatchField;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            faPatchField,
            patch,
            (
                const faPatch& p,
                const DimensionedField<Type, areaMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            faPatchField,
            patchMapper,
            (
                const faPatchField<Type>& ptf,
                const faPatch& p,
                const DimensionedField<Type, areaMesh>& iF,
                const faPatchFieldMapper& m
            ),
            (dynamic_cast<const faPatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            faPatchField,
            dictionary,
            (
                const faPatch& p,
                const DimensionedField<Type, areaMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        faPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch and internal field and patch field
        faPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        faPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given faPatchField onto a new patch
        faPatchField
        (
            const faPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        faPatchField(const faPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<faPatchField<Type> > clone() const
        {
            return tmp<faPatchField<Type> >(new faPatchField<Type>(*this));
        }

        //- Construct as copy setting internal field reference
        faPatchField
        (
            const faPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<Type> > clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return tmp<faPatchField<Type> >(new faPatchField<Type>(*this, iF));
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<faPatchField<Type> > New
        (
            const word&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given faPatchField mapped onto a new patch
        static tmp<faPatchField<Type> > New
        (
            const faPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<faPatchField<Type> > New
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedFaPatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static tmp<faPatchField<Type> > NewCalculatedType
        (
            const faPatchField<Type2>&
        );


    // Destructor

        virtual ~faPatchField<Type>()
        {}


    // Member functions

        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return patch
            const faPatch& patch() const
            {
                return patch_;
            }

            //- Return dimensioned internal field reference
            const DimensionedField<Type, areaMesh>&
            dimensionedInternalField() const
            {
                return internalField_;
            }

            //- Return internal field reference
            const Field<Type>& internalField() const
            {
                return internalField_;
            }

            //- Return the type of the calculated for of faPatchField
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

            //- Return true if the boundary condition has already been updated
            bool updated() const
            {
                return updated_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const faPatchFieldMapper&
            );

            //- Reverse map the given faPatchField onto this faPatchField
            virtual void rmap
            (
                const faPatchField<Type>&,
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
                const Pstream::commsTypes commsType=Pstream::blocking
            )
            {}

            //- Evaluate the patch field, sets Updated to false
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
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


        //- Write
        virtual void write(Ostream&) const;


        // Check

            //- Check faPatchField<Type> against given faPatchField<Type>
            void check(const faPatchField<Type>&) const;


    // Member operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const faPatchField<Type>&);
        virtual void operator+=(const faPatchField<Type>&);
        virtual void operator-=(const faPatchField<Type>&);
        virtual void operator*=(const faPatchField<scalar>&);
        virtual void operator/=(const faPatchField<scalar>&);

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

        virtual void operator==(const faPatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);


    // Ostream operator

        friend Ostream& operator<< <Type>(Ostream&, const faPatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faPatchField.C"
#   include "calculatedFaPatchField.H"
#endif


#define makeFaPatchTypeFieldTypeName(typePatchTypeField)                   \
                                                                           \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);

#define makeFaPatchFieldsTypeName(typePatchField)                          \
                                                                           \
makeFaPatchTypeFieldTypeName(typePatchField##FaPatchScalarField);          \
makeFaPatchTypeFieldTypeName(typePatchField##FaPatchVectorField);          \
makeFaPatchTypeFieldTypeName(typePatchField##FaPatchSphericalTensorField); \
makeFaPatchTypeFieldTypeName(typePatchField##FaPatchSymmTensorField);      \
makeFaPatchTypeFieldTypeName(typePatchField##FaPatchTensorField);

#define makeFaPatchTypeField(PatchTypeField, typePatchTypeField)           \
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


#define makeFaPatchFields(type)                                            \
                                                                           \
makeFaPatchTypeField(faPatchScalarField, type##FaPatchScalarField);        \
makeFaPatchTypeField(faPatchVectorField, type##FaPatchVectorField);        \
makeFaPatchTypeField                                                       \
(                                                                          \
    faPatchSphericalTensorField,                                           \
    type##FaPatchSphericalTensorField                                      \
);                                                                         \
makeFaPatchTypeField(faPatchSymmTensorField, type##FaPatchSymmTensorField);\
makeFaPatchTypeField(faPatchTensorField, type##FaPatchTensorField);


#define makeFaPatchTypeFieldTypedefs(type)                                   \
                                                                             \
typedef type##FaPatchField<scalar> type##FaPatchScalarField;                 \
typedef type##FaPatchField<vector> type##FaPatchVectorField;                 \
typedef type##FaPatchField<sphericalTensor>                                  \
     type##FaPatchSphericalTensorField;                                      \
typedef type##FaPatchField<symmTensor> type##FaPatchSymmTensorField;         \
typedef type##FaPatchField<tensor> type##FaPatchTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
