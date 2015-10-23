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
    chtRegionCoupleBase

Description
    Region couple patch field taking radiation into account

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved

SourceFiles
    chtRegionCoupleBase.C

\*---------------------------------------------------------------------------*/

#ifndef chtRegionCoupleBase_H
#define chtRegionCoupleBase_H

#include "regionCouplingFvPatchFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class chtRcTemperatureFvPatchScalarField;

/*---------------------------------------------------------------------------*\
                  Class chtRegionCoupleBase Declaration
\*---------------------------------------------------------------------------*/

class chtRegionCoupleBase
:
    public regionCouplingFvPatchScalarField
{
    // Private data


public:

    //- Runtime type information
    TypeName("chtRegionCoupleBase");


    // Constructors

        //- Construct from patch and internal field
        chtRegionCoupleBase
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        chtRegionCoupleBase
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given chtRegionCoupleBase onto a new patch
        chtRegionCoupleBase
        (
            const chtRegionCoupleBase&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        chtRegionCoupleBase
        (
            const chtRegionCoupleBase&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new chtRegionCoupleBase(*this)
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
                new chtRegionCoupleBase(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return false to disable interpolation on this patch
            virtual bool coupled() const
            {
                return false;
            }

            //- Return shadow patch field
            const chtRegionCoupleBase& shadowPatchField() const;

            //- Return original patch field
            const scalarField& originalPatchField() const
            {
                return regionCouplingFvPatchScalarField::originalPatchField();
            }


        // Evaluation functions

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

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

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
