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
    engineMassFlowRateInletVelocityFvPatchVectorField

Description
    Describes an massflow normal vector boundary condition by its magnitude
    as an integral over its area.
    The current density is used to correct the velocity.
    Note: The value is positive inwards

    Example of the BC specification:

        inlet
        {
            type            massFlowRateInletVelocity;
            massFlowRate    0.2;     // Mass flow rate [kg/s]
        }

NOTE!
 * may not work correctly for transonic inlets!
 * strange behaviour with potentialFoam since the U equation is not solved

SourceFiles
    engineMassFlowRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef engineMassFlowRateInletVelocityFvPatchVectorField_H
#define engineMassFlowRateInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "engineTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
          Class massFlowRateInletVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/
class graph;

class engineMassFlowRateInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalize the mass flux
        word rhoName_;

        fileName timeDataFileName_;
        autoPtr<graph> timeDataPtr_;

        const engineTime& engineDB_;

        void checkTable();


public:

   //- Runtime type information
   TypeName("engineMassFlowRateInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        engineMassFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        engineMassFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  engineMassFlowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        engineMassFlowRateInletVelocityFvPatchVectorField
        (
            const engineMassFlowRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        engineMassFlowRateInletVelocityFvPatchVectorField
        (
            const engineMassFlowRateInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new engineMassFlowRateInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        engineMassFlowRateInletVelocityFvPatchVectorField
        (
            const engineMassFlowRateInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new engineMassFlowRateInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the fluctuation scale
            const graph& timeData() const
            {
                return timeDataPtr_;
            }

/*
            //- Return reference to the fluctuation scale to allow adjustment
            graph& timeData()
            {
                return timeDataPtr_;
            }
*/
            const engineTime& engineDB() const
            {
                return engineDB_;
            }


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
