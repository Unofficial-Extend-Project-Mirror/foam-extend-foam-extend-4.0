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
    solidWallMixedTemperatureCoupledFvPatchScalarField

Description
    Mixed boundary condition for temperature, to be used by the
    conjugate heat transfer solver.
    If my temperature is T1, neighbour is T2:

    T1 > T2: my side becomes fixedValue T2 bc, other side becomes fixedGradient.
    

    Example usage:
        myInterfacePatchName
        {
            type                solidWallMixedTemperatureCoupled;
            neighbourRegionName fluid;
            neighbourPatchName  fluidSolidInterface;
            neighbourFieldName  T;
            K                   K;
            value               uniform 300;
        }

SourceFiles
    solidWallMixedTemperatureCoupledFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef solidWallMixedTemperatureCoupledFvPatchScalarField_H
#define solidWallMixedTemperatureCoupledFvPatchScalarField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "coupleManager.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class solidWallMixedTemperatureCoupledFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class solidWallMixedTemperatureCoupledFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Couple manager object
        coupleManager coupleManager_;

        //- Name of thermal conductivity field
        word KName_;

        bool fixesValue_;

public:

    //- Runtime type information
    TypeName("solidWallMixedTemperatureCoupled");


    // Constructors

        //- Construct from patch and internal field
        solidWallMixedTemperatureCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidWallMixedTemperatureCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidWallMixedTemperatureCoupledFvPatchScalarField
        //  onto a new patch
        solidWallMixedTemperatureCoupledFvPatchScalarField
        (
            const solidWallMixedTemperatureCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new solidWallMixedTemperatureCoupledFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidWallMixedTemperatureCoupledFvPatchScalarField
        (
            const solidWallMixedTemperatureCoupledFvPatchScalarField&,
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
                new solidWallMixedTemperatureCoupledFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Get corresponding K field
        const fvPatchScalarField& K() const;

        //- Return true if this patch field fixes a value.
        //  Needed to check if a level has to be specified while solving
        //  Poissons equations.
        virtual bool fixesValue() const
        {
            return fixesValue_;
        }

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::blocking
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
