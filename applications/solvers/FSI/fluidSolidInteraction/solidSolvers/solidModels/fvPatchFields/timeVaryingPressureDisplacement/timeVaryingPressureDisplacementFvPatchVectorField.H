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
    timeVaryingPressureDisplacementFvPatchVectorField

Description
    Fixed traction boundary condition for the standard linear elastic, fixed
    coefficient displacement equation (stressedFoam).

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    timeVaryingPressureDisplacementFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingPressureDisplacementFvPatchVectorField_H
#define timeVaryingPressureDisplacementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "tractionDisplacementFvPatchVectorField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class timeVaryingPressureDisplacementFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingPressureDisplacementFvPatchVectorField
:
    public tractionDisplacementFvPatchVectorField
{

    // Private Data

        //- The time series being used, including the bounding treatment
        interpolationTable<scalar> timeSeries_;

public:

    //- Runtime type information
    TypeName("timeVaryingPressureDisplacement");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingPressureDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingPressureDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  timeVaryingPressureDisplacementFvPatchVectorField onto a new patch
        timeVaryingPressureDisplacementFvPatchVectorField
        (
            const timeVaryingPressureDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingPressureDisplacementFvPatchVectorField
        (
            const timeVaryingPressureDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingPressureDisplacementFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingPressureDisplacementFvPatchVectorField
        (
            const timeVaryingPressureDisplacementFvPatchVectorField&,
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
                new timeVaryingPressureDisplacementFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return the time series used
            const interpolationTable<scalar>& timeSeries() const
            {
                return timeSeries_;
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
