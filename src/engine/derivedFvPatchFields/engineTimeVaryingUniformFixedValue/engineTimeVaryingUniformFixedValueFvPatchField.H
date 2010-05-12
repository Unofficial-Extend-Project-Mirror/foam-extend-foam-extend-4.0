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
    engineTimeVaryingUniformFixedValueFvPatchField

Description

SourceFiles
    engineTimeVaryingUniformFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef engineTimeVaryingUniformFixedValueFvPatchField_H
#define engineTimeVaryingUniformFixedValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "engineTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class graph;

/*---------------------------------------------------------------------------*\
         Class engineTimeVaryingUniformFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class engineTimeVaryingUniformFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        fileName timeDataFileName_;
        autoPtr<graph> timeDataPtr_;

        const engineTime& engineDB_;

        void checkTable();


public:

    //- Runtime type information
    TypeName("engineTimeVaryingUniformFixedValue");


    // Constructors

        //- Construct from patch and internal field
        engineTimeVaryingUniformFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        engineTimeVaryingUniformFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given engineTimeVaryingUniformFixedValueFvPatchField
        //  onto a new patch
        engineTimeVaryingUniformFixedValueFvPatchField
        (
            const engineTimeVaryingUniformFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        engineTimeVaryingUniformFixedValueFvPatchField
        (
            const engineTimeVaryingUniformFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new engineTimeVaryingUniformFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        engineTimeVaryingUniformFixedValueFvPatchField
        (
            const engineTimeVaryingUniformFixedValueFvPatchField<Type>&,
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
                new engineTimeVaryingUniformFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the fluctuation scale
            const graph& timeData() const
            {
                return timeDataPtr_;
            }

            //- Return reference to the fluctuation scale to allow adjustment
            graph& timeData()
            {
                return timeDataPtr_;
            }

            const engineTime& engineDB() const
            {
                return engineDB_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};

template<>
void engineTimeVaryingUniformFixedValueFvPatchField<scalar>::updateCoeffs();

template<>
void engineTimeVaryingUniformFixedValueFvPatchField<vector>::updateCoeffs();

template<>
void engineTimeVaryingUniformFixedValueFvPatchField<sphericalTensor>::updateCoeffs();

template<>
void engineTimeVaryingUniformFixedValueFvPatchField<symmTensor>::updateCoeffs();

template<>
void engineTimeVaryingUniformFixedValueFvPatchField<tensor>::updateCoeffs();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "engineTimeVaryingUniformFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
