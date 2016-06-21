/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    velocityTractionDisplacementFvPatchVectorField

Description
    Fixed traction boundary condition for the linear/nonlinear elastic, fixed
    coefficient displacement equation.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

SourceFiles
    velocityTractionDisplacementFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef velocityTractionDisplacementFvPatchVectorField_H
#define velocityTractionDisplacementFvPatchVectorField_H

// #include "fvPatchFields.H"
#include "robinFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class velocityTractionDisplacementFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class velocityTractionDisplacementFvPatchVectorField
:
    public robinFvPatchVectorField
{

    // Private Data

        //- Traction
        vectorField traction_;

        //- Pressure
        scalarField pressure_;

        //- Velocity
        vectorField velocity_;

        //- Current normal
        vectorField normal_;

        //- Second order accuracy correction
        Switch secondOrder_;

        //- Limiter coefficient for snGrad correction
        scalar limitCoeff_;

public:

    //- Runtime type information
    TypeName("velocityTractionDisplacement");


    // Constructors

        //- Construct from patch and internal field
        velocityTractionDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        velocityTractionDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  velocityTractionDisplacementFvPatchVectorField onto a new patch
        velocityTractionDisplacementFvPatchVectorField
        (
            const velocityTractionDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        velocityTractionDisplacementFvPatchVectorField
        (
            const velocityTractionDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new velocityTractionDisplacementFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        velocityTractionDisplacementFvPatchVectorField
        (
            const velocityTractionDisplacementFvPatchVectorField&,
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
                new velocityTractionDisplacementFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            virtual const vectorField& traction() const
            {
                return traction_;
            }

            virtual vectorField& traction()
            {
                return traction_;
            }

            virtual const scalarField& pressure() const
            {
                return pressure_;
            }

            virtual  scalarField& pressure()
            {
                return pressure_;
            }

            virtual const vectorField& velocit() const
            {
                return velocity_;
            }

            virtual vectorField& velocity()
            {
                return velocity_;
            }

            virtual const vectorField& normal() const
            {
                return normal_;
            }

            virtual vectorField& normal()
            {
                return normal_;
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
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

//         //- Evaluate the patch field
//         virtual void evaluate
//         (
//             const Pstream::commsTypes commsType=Pstream::blocking
//         );

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

//         virtual void operator=(const fvPatchField<vector>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
