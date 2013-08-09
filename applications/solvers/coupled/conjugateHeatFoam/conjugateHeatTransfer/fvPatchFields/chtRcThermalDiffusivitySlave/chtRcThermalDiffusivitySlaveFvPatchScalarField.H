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
    chtRcThermalDiffusivitySlaveFvPatchScalarField

Description
    Region couple patch field taking radiation into account

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved

SourceFiles
    chtRcThermalDiffusivitySlaveFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef chtRcThermalDiffusivitySlaveFvPatchScalarField_H
#define chtRcThermalDiffusivitySlaveFvPatchScalarField_H

#include "chtRegionCoupleBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class chtThermalDiffusivityRegionCoupleFvPatchScalarField;
class chtTemperatureRegionCoupleFvPatchScalarField;

/*---------------------------------------------------------------------------*\
                  Class chtRcThermalDiffusivitySlaveFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class chtRcThermalDiffusivitySlaveFvPatchScalarField
:
    public chtRegionCoupleBase
{
    // Private data

public:

    //- Runtime type information
    TypeName("chtRcThermalDiffusivitySlave");


    // Constructors

        //- Construct from patch and internal field
        chtRcThermalDiffusivitySlaveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        chtRcThermalDiffusivitySlaveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given chtRcThermalDiffusivitySlaveFvPatchScalarField onto a new patch
        chtRcThermalDiffusivitySlaveFvPatchScalarField
        (
            const chtRcThermalDiffusivitySlaveFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        chtRcThermalDiffusivitySlaveFvPatchScalarField
        (
            const chtRcThermalDiffusivitySlaveFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new chtRcThermalDiffusivitySlaveFvPatchScalarField(*this)
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
                new chtRcThermalDiffusivitySlaveFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return shadow patch field
            virtual const chtRegionCoupleBase& shadowPatchField() const;


        // Evaluation functions

            //- Evaluate the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate(const Pstream::commsTypes commsType);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
