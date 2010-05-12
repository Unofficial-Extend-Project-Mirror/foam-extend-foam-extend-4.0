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
    Foam::incompressible::RASModels::nutStandardRoughWallFunctionFvPatchScalarField

Description
    Wall function boundary condition for rough walls

SourceFiles
    nutStandardWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutStandardRoughWallFunctionFvPatchScalarField_H
#define nutStandardRoughWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                    Class nutWallFunctionFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nutStandardRoughWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;

        scalar roughnessHeight_;
        scalar roughnessConstant_;
        scalar roughnessFudgeFactor_;


public:

    //- Runtime type information
    TypeName("nutStandardRoughWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutStandardRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutStandardRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutStandardRoughWallFunctionFvPatchScalarField
        //  onto a new patch
        nutStandardRoughWallFunctionFvPatchScalarField
        (
            const nutStandardRoughWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutStandardRoughWallFunctionFvPatchScalarField
        (
            const nutStandardRoughWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutStandardRoughWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutStandardRoughWallFunctionFvPatchScalarField
        (
            const nutStandardRoughWallFunctionFvPatchScalarField&,
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
                new nutStandardRoughWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of the velocity field
            const word& UName() const
            {
                return UName_;
            }


            //- Return the fluctuation scale
            const scalar& roughnessHeight() const
            {
                return roughnessHeight_;
            }

            //- Return reference to the fluctuation scale to allow adjustment
            scalar& roughnessHeight()
            {
                return roughnessHeight_;
            }


            //- Return the fluctuation scale
            const scalar& roughnessConstant() const
            {
                return roughnessConstant_;
            }

            //- Return reference to the fluctuation scale to allow adjustment
            scalar& roughnessConstant()
            {
                return roughnessConstant_;
            }

            //- Return the fluctuation scale
            const scalar& roughnessFudgeFactor() const
            {
                return roughnessFudgeFactor_;
            }

            //- Return reference to the fluctuation scale to allow adjustment
            scalar& roughnessFudgeFactor()
            {
                return roughnessFudgeFactor_;
            }


        // Evaluation functions

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );


        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
