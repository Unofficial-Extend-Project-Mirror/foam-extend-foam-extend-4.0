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
    Foam::wideBandSpecularRadiationMixedFvPatchScalarField

Description
    Foam::wideBandSpecularRadiationMixedFvPatchScalarField

SourceFiles
    wideBandSpecularRadiationMixedFvPatchScalarField.C

Author
    Dominik Christ, Wikki Ltd. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef wideBandSpecularRadiationMixedFvPatchScalarField_H
#define wideBandSpecularRadiationMixedFvPatchScalarField_H

#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
            Class wideBandSpecularRadiationMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wideBandSpecularRadiationMixedFvPatchScalarField
:
    public mixedFvPatchField<scalar>
{

    //- For each face, ray index of ray in direction of reflection
    mutable labelListList *receivedRayIDPtr_;

    //- Calculate IDs of rays in direction of reflection
    void calcReceivedRayIDs() const;

public:

    //- Runtime type information
    TypeName("wideBandSpecularRadiation");


    // Constructors

        //- Construct from patch and internal field
        wideBandSpecularRadiationMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wideBandSpecularRadiationMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given wideBandSpecularRadiationMixedFvPatchScalarField onto a new patch
        wideBandSpecularRadiationMixedFvPatchScalarField
        (
            const wideBandSpecularRadiationMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wideBandSpecularRadiationMixedFvPatchScalarField
        (
            const wideBandSpecularRadiationMixedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new wideBandSpecularRadiationMixedFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wideBandSpecularRadiationMixedFvPatchScalarField
        (
            const wideBandSpecularRadiationMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new wideBandSpecularRadiationMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return ray IDs of reflected rays
            const labelListList& receivedRayIDs() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchScalarField& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
