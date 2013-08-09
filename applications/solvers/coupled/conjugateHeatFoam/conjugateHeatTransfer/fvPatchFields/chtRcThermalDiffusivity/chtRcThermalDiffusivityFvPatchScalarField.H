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
    chtRcThermalDiffusivityFvPatchScalarField

Description
    Region couple patch field taking radiation into account

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved

SourceFiles
    chtRcThermalDiffusivityFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef chtRcThermalDiffusivityFvPatchScalarField_H
#define chtRcThermalDiffusivityFvPatchScalarField_H

#include "chtRegionCoupleBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class chtRcThermalDiffusivityFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class chtRcThermalDiffusivityFvPatchScalarField
:
    public chtRegionCoupleBase
{
    // Private data


public:

    //- Runtime type information
    TypeName("chtRcThermalDiffusivity");


    // Constructors

        //- Construct from patch and internal field
        chtRcThermalDiffusivityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        chtRcThermalDiffusivityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given chtRcThermalDiffusivityFvPatchScalarField onto a new patch
        chtRcThermalDiffusivityFvPatchScalarField
        (
            const chtRcThermalDiffusivityFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        chtRcThermalDiffusivityFvPatchScalarField
        (
            const chtRcThermalDiffusivityFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new chtRcThermalDiffusivityFvPatchScalarField(*this)
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
                new chtRcThermalDiffusivityFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access


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
