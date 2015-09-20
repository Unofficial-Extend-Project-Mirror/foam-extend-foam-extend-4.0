/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Foam::greyDiffusiveRadiationMixedFvPatchScalarField

Description
    Radiation boundary will diffusive reflection and black body emission.
    See Modest "Radiative heat transfer", Capter 16.6 "The finite volume method"
    Eq.16.63

SourceFiles
    greyDiffusiveRadiationMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef greyDiffusiveRadiationMixedFvPatchScalarField_H
#define greyDiffusiveRadiationMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
/*---------------------------------------------------------------------------*\
        Class greyDiffusiveRadiationMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class greyDiffusiveRadiationMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of temperature field
        word TName_;

        //- Emissivity
        scalar emissivity_;

        //- Sum of cosine factors of all outgoing rays (lazy evaluated)
        mutable scalarField *sumOutgoingAnglesPtr_;


    // Private member functions

        //- Calculate sum of cosine factors of all outgoing rays
        void calcSumOutgoingAngles() const;


public:

    //- Runtime type information
    TypeName("greyDiffusiveRadiation");


    // Constructors

        //- Construct from patch and internal field
        greyDiffusiveRadiationMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        greyDiffusiveRadiationMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        //  greyDiffusiveRadiationMixedFvPatchScalarField onto a new patch
        greyDiffusiveRadiationMixedFvPatchScalarField
        (
            const greyDiffusiveRadiationMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        greyDiffusiveRadiationMixedFvPatchScalarField
        (
            const greyDiffusiveRadiationMixedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new greyDiffusiveRadiationMixedFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        greyDiffusiveRadiationMixedFvPatchScalarField
        (
            const greyDiffusiveRadiationMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new greyDiffusiveRadiationMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the temperature field name
            const word& TName() const
            {
                return TName_;
            }

            //- Return reference to the temperature field name to allow
            //  adjustment
            word& TName()
            {
                return TName_;
            }

            //- Return the emissivity
            const scalar& emissivity() const
            {
                return emissivity_;
            }

            //- Return reference to the emissivity to allow adjustment
            scalar& emissivity()
            {
                return emissivity_;
            }

            //- Return sum of cosine factors of all outgoing rays
            const scalarField& sumOutgoingAngles() const;


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
