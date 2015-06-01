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
    Foam::alphaContactAngleFvPatchScalarField

Description
    Abstract base class for alphaContactAngle boundary conditions.

    Derived classes must implement the theta() fuction which returns the
    wall contact angle field.

SourceFiles
    alphaContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphaContactAngleFvPatchScalarField_H
#define alphaContactAngleFvPatchScalarField_H

#include "zeroGradientFvPatchFields.H"
#include "fvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class alphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class alphaContactAngleFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("alphaContactAngle");


    // Constructors

        //- Construct from patch and internal field
        alphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given alphaContactAngleFvPatchScalarField
        //  onto a new patch
        alphaContactAngleFvPatchScalarField
        (
            const alphaContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphaContactAngleFvPatchScalarField
        (
            const alphaContactAngleFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        alphaContactAngleFvPatchScalarField
        (
            const alphaContactAngleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member functions

        //- Return the contact angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
