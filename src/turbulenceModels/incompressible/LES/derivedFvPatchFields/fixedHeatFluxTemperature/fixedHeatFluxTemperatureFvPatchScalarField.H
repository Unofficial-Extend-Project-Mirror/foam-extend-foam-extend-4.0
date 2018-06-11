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
    Foam::fixedHeatFluxTemperatureFvPatchScalarField

Description

    Fixed temperature flux boundary condition. Assumes incompressible turbulent
    flow. The normal temperature gradient is calculated as:

    snGrad(T) = q/(kappaEff*rhoRef*c)

    where:
        - q is the prescribed heat flux [W/m^2],
        - kappaEff is the effective diffusivity field: nu/Pr + nu_t/Pr_t,
        - Pr and Prt are laminar and turbulent Prandtl number read from
          transportProperties dictionary
        - rhoRef is the reference density
        - c is the specific heat capacity

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved

SourceFiles
    fixedHeatFluxTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedHeatFluxTemperatureFvPatchScalarField_H
#define fixedHeatFluxTemperatureFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class fixedHeatFluxTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedHeatFluxTemperatureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Fixed heat flux
        scalarField heatFlux_;

        //- Laminar Prandtl number
        scalar Pr_;

        //- Turbulent Prandtl number
        scalar Prt_;

        //- Reference density
        scalar rhoRef_;

        //- Specific heat capacity of the fluid
        scalar c_;


public:

    //- Runtime type information
    TypeName("fixedHeatFluxTemperature");


    // Constructors

        //- Construct from patch and internal field
        fixedHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedHeatFluxTemperatureFvPatchScalarField onto a new patch
        fixedHeatFluxTemperatureFvPatchScalarField
        (
            const fixedHeatFluxTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedHeatFluxTemperatureFvPatchScalarField
        (
            const fixedHeatFluxTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedHeatFluxTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedHeatFluxTemperatureFvPatchScalarField
        (
            const fixedHeatFluxTemperatureFvPatchScalarField&,
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
                new fixedHeatFluxTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
