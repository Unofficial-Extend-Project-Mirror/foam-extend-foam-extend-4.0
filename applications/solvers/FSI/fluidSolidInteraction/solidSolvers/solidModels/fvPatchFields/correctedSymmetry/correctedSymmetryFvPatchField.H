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
    Foam::correctedSymmetryFvPatchField

Description
    Foam::correctedSymmetryFvPatchField

SourceFiles
    correctedSymmetryFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef correctedSymmetryFvPatchField_H
#define correctedSymmetryFvPatchField_H

#include "symmetryFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class correctedSymmetryFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class correctedSymmetryFvPatchField
:
    public symmetryFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("correctedSymmetry");
//     TypeName(symmetryFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        correctedSymmetryFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        correctedSymmetryFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  correctedSymmetryFvPatchField onto a new patch
        correctedSymmetryFvPatchField
        (
            const correctedSymmetryFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        correctedSymmetryFvPatchField
        (
            const correctedSymmetryFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new correctedSymmetryFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        correctedSymmetryFvPatchField
        (
            const correctedSymmetryFvPatchField<Type>&,
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
                new correctedSymmetryFvPatchField<Type>(*this, iF)
            );
        }

    // Member functions

        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type> > snGrad() const
            {
                return symmetryFvPatchField<Type>::snGrad();
            }

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            )
            {
                symmetryFvPatchField<Type>::evaluate(commsType);
            }
};

// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

// template<>
// tmp<scalarField > correctedSymmetryFvPatchField<scalar>::snGrad() const;

template<>
void correctedSymmetryFvPatchField<scalar>::evaluate
(
    const Pstream::commsTypes commsType
);

template<>
tmp<vectorField > correctedSymmetryFvPatchField<vector>::snGrad() const;

template<>
void correctedSymmetryFvPatchField<vector>::evaluate
(
    const Pstream::commsTypes commsType
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "correctedSymmetryFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
