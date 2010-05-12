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
    solidWallHeatFluxCoupledFvPatchScalarField

Description
    Fixed value boundary condition for temperature, to be used by the
    conjugate heat transfer solver.

    Example usage:
        myInterfacePatchName
        {
            type                solidWallTemperatureCoupled;
            neighbourRegionName fluid;
            neighbourPatchName  fluidSolidInterface;
            neighbourFieldName  T;
            K                   K;
            value               uniform 300;
        }

SourceFiles
    solidWallTemperatureCoupledFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef solidWallTemperatureCoupledFvPatchScalarField_H
#define solidWallTemperatureCoupledFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "coupleManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class solidWallTemperatureCoupledFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class solidWallTemperatureCoupledFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Couple manager object
        coupleManager coupleManager_;

        //- Name of thermal conductivity field
        word KName_;


public:

    //- Runtime type information
    TypeName("solidWallTemperatureCoupled");


    // Constructors

        //- Construct from patch and internal field
        solidWallTemperatureCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidWallTemperatureCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidWallTemperatureCoupledFvPatchScalarField
        //  onto a new patch
        solidWallTemperatureCoupledFvPatchScalarField
        (
            const solidWallTemperatureCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new solidWallTemperatureCoupledFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidWallTemperatureCoupledFvPatchScalarField
        (
            const solidWallTemperatureCoupledFvPatchScalarField&,
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
                new solidWallTemperatureCoupledFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- (intensive) flux
        tmp<scalarField> flux() const;

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
