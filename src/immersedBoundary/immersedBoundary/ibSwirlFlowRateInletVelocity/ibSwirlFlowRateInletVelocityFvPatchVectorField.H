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
    Foam::ibSwirlFlowRateInletVelocityFvPatchVectorField

Description
    Describes a volumetric/mass flow normal vector boundary condition by its
    magnitude as an integral over its area with a swirl component determined
    by the RPM

    The basis of the patch (volumetric or mass) is determined by the
    dimensions of the flux, phi.
    The current density is used to correct the velocity when applying the
    mass basis.

    The boundary condition is sensitised to work only on live IB cells

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            ibSwirlFlowRateInletVelocity;
        flowRate        0.2;        // Volumetric/mass flow rate [m3/s or kg/s]
        rpm             100;
    }
    \endverbatim

Note
    - The value is positive inwards

SourceFiles
    ibSwirlFlowRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef ibSwirlFlowRateInletVelocityFvPatchVectorField_H
#define ibSwirlFlowRateInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
       Class ibSwirlFlowRateInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class ibSwirlFlowRateInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet integral flow rate
        const scalar flowRate_;

        //- Name of the flux transporting the field
        const word phiName_;

        //- Name of the density field used to normalize the mass flux
        const word rhoName_;

        //- RPM
        const scalar rpm_;

        // Derived parameters

            //- Total area of inlet patch
            scalar totArea_;

            //- Average centre
            vector avgCenter_;

            //- Average normal
            vector avgNormal_;


    // Private member functions

        //- Calculate geometrical characteristics of the patch
        void calcGeom();


public:

   //- Runtime type information
   TypeName("ibSwirlFlowRateInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        ibSwirlFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ibSwirlFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        ibSwirlFlowRateInletVelocityFvPatchVectorField
        (
            const ibSwirlFlowRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ibSwirlFlowRateInletVelocityFvPatchVectorField
        (
            const ibSwirlFlowRateInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new ibSwirlFlowRateInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        ibSwirlFlowRateInletVelocityFvPatchVectorField
        (
            const ibSwirlFlowRateInletVelocityFvPatchVectorField&,
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
                new ibSwirlFlowRateInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the flux
            scalar flowRate() const
            {
                return flowRate_;
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
