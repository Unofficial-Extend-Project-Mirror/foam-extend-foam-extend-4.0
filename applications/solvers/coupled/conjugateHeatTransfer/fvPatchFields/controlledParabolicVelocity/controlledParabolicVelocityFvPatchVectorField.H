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
    controlledParabolicVelocityFvPatchVectorField

Description
    Boundary condition specifies a velocity inlet profile
    (fixed value), given the volumetric flow rate
    flow direction n and direction of the parabolic coordinate y

SourceFiles
    controlledParabolicVelocityFvPatchVectorField.C

Author
    Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef controlledParabolicVelocityFvPatchVectorField_H
#define controlledParabolicVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class parabolicVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class controlledParabolicVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Volumetric flow rate
        scalar Umean_;

        //- Flow direction
        vector n_;

        //- 1D coordinate direction
        vector y_;

        //- target value
        scalar target_;

        //- field to observe
        word obsFieldName_;

        //- patch ID to observe
        word obsPatchName_;

        //- patch ID to observe
        label obsPatchID_;

        //- proportionality factor
        scalar gain_;

        //- Current time index
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("controlledParabolicVelocity");


    // Constructors

        //- Construct from patch and internal field
        controlledParabolicVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        controlledParabolicVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given controlledParabolicVelocityFvPatchVectorField
        //  onto a new patch
        controlledParabolicVelocityFvPatchVectorField
        (
            const controlledParabolicVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new controlledParabolicVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        controlledParabolicVelocityFvPatchVectorField
        (
            const controlledParabolicVelocityFvPatchVectorField&,
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
                new controlledParabolicVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
