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
    regionCouplingResistanceFvPatchField

Description
    Region couple patch field jump resistance on the surface taking
    radiation into account

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved

SourceFiles
    regionCouplingResistanceFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCouplingResistanceFvPatchField_H
#define regionCouplingResistanceFvPatchField_H

#include "chtRegionCoupleBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class regionCouplingResistanceFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class chtRcThermalDiffusivityResistanceFvPatchScalarField
:
    public chtRegionCoupleBase
{
    // Private data

        //- Conductivity
        scalarField conductivity_;


public:

    //- Runtime type information
    TypeName("chtRcThermalDiffusivityResistance");


    // Constructors

        //- Construct from patch and internal field
        chtRcThermalDiffusivityResistanceFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        chtRcThermalDiffusivityResistanceFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCouplingResistanceFvPatchField
        //  onto a new patch
        chtRcThermalDiffusivityResistanceFvPatchScalarField
        (
            const chtRcThermalDiffusivityResistanceFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        chtRcThermalDiffusivityResistanceFvPatchScalarField
        (
            const chtRcThermalDiffusivityResistanceFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new chtRcThermalDiffusivityResistanceFvPatchScalarField(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new chtRcThermalDiffusivityResistanceFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return conductivity
            const scalarField& conductivity() const
            {
                return conductivity_;
            }

            //- Return Access to conductivity
            scalarField& conductivity()
            {
                return conductivity_;
            }


        // Evaluation functions

            //- Evaluate the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate(const Pstream::commsTypes commsType);

            //- Calculate the thermal diffusivity
            virtual void calcThermalDiffusivity
            (
                chtRegionCoupleBase& owner,
                const chtRegionCoupleBase& neighbour
            ) const;

            //- Calculate the temperature
            virtual void calcTemperature
            (
                chtRcTemperatureFvPatchScalarField& owner,
                const chtRcTemperatureFvPatchScalarField& neighbour,
                const chtRegionCoupleBase& ownerK
            ) const;


        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
