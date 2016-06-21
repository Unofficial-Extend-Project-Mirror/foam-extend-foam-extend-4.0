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
    solidTractionFvPatchVectorField

Description
    Fixed traction boundary condition for solid foam solvers.
    Suitable for for total strain or incremental strain total or
    updated Lagrangian approaches.
    Contains thermal-elastic, plastic, non-linear terms if required.
    And non-orthogonal correction.

Author
   Philip Cardiff

SourceFiles
    solidTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidTractionFvPatchVectorField_H
#define solidTractionFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "nonLinearGeometry.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class solidTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidTractionFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{

    // Private Data

        //- Name of the displacement field
        word fieldName_;

        //- Traction
        vectorField traction_;

        //- Pressure
        scalarField pressure_;

        //- Type of geometrical non-linearity
       nonLinearGeometry::nonLinearType nonLinear_;

        //- if it is an orthropic solver
        Switch orthotropic_;


public:

    //- Runtime type information
    TypeName("solidTraction");


    // Constructors

        //- Construct from patch and internal field
        solidTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  solidTractionFvPatchVectorField onto a new patch
        solidTractionFvPatchVectorField
        (
            const solidTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidTractionFvPatchVectorField
        (
            const solidTractionFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new solidTractionFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidTractionFvPatchVectorField
        (
            const solidTractionFvPatchVectorField&,
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
                new solidTractionFvPatchVectorField(*this, iF)
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

            virtual scalarField& pressure()
            {
                return pressure_;
            }

            virtual word& fieldName()
            {
                return fieldName_;
            }

            const nonLinearGeometry::nonLinearType& nonLinear() const
            {
                return nonLinear_;
            }

            nonLinearGeometry::nonLinearType& nonLinear()
            {
                return nonLinear_;
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
