/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | 
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
    Foam::realFluidWallHeatTransferFvPatchScalarField

Description
    Enthalpy boundary conditions for wall heat transfer

SourceFiles
    realFluidWallHeatTransferFvPatchScalarField.C

Modified by
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#ifndef realFluidWallHeatTransferFvPatchScalarField_H
#define realFluidWallHeatTransferFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class wallHeatTransferFvPatch Declaration
\*---------------------------------------------------------------------------*/

class realFluidWallHeatTransferFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Tinf
        scalarField Tinf_;

        //- alphaWall
        scalarField alphaWall_;


public:

    //- Runtime type information
    TypeName("realFluidWallHeatTransfer");


    // Constructors

        //- Construct from patch and internal field
        realFluidWallHeatTransferFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        realFluidWallHeatTransferFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given realFluidWallHeatTransferFvPatchScalarField
        // onto a new patch
        realFluidWallHeatTransferFvPatchScalarField
        (
            const realFluidWallHeatTransferFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        realFluidWallHeatTransferFvPatchScalarField
        (
            const realFluidWallHeatTransferFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new realFluidWallHeatTransferFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        realFluidWallHeatTransferFvPatchScalarField
        (
            const realFluidWallHeatTransferFvPatchScalarField&,
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
                new realFluidWallHeatTransferFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return Tinf
            const scalarField& Tinf() const
            {
                return Tinf_;
            }

            //- Return reference to Tinf to allow adjustment
            scalarField& Tinf()
            {
                return Tinf_;
            }

            //- Return alphaWall
            const scalarField& alphaWall() const
            {
                return alphaWall_;
            }

            //- Return reference to alphaWall to allow adjustment
            scalarField& alphaWall()
            {
                return alphaWall_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

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
