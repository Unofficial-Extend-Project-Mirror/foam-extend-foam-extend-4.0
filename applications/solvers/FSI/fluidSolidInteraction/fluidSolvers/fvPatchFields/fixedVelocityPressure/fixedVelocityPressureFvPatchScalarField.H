/*---------------------------------------------------------------------------* \
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
    Foam::fixedVelocityPressureFvPatchScalarField

Description
    Set the pressure gradient boundary condition appropriately for boundary where
    velocity is fixed.

SourceFiles
    fixedVelocityPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedVelocityPressureFvPatchScalarFields_H
#define fixedVelocityPressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fixedVelocityPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedVelocityPressureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

public:

    //- Runtime type information
    TypeName("fixedVelocityPressure");


    // Constructors

        //- Construct from patch and internal field
        fixedVelocityPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedVelocityPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedVelocityPressureFvPatchScalarField onto a new patch
        fixedVelocityPressureFvPatchScalarField
        (
            const fixedVelocityPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedVelocityPressureFvPatchScalarField
        (
            const fixedVelocityPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedVelocityPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedVelocityPressureFvPatchScalarField
        (
            const fixedVelocityPressureFvPatchScalarField&,
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
                new fixedVelocityPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //-
        virtual void updateCoeffs();

//         //- Evaluate the patch field
//         virtual void evaluate
//         (
//             const Pstream::commsTypes commsType=Pstream::blocking
//         );


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
