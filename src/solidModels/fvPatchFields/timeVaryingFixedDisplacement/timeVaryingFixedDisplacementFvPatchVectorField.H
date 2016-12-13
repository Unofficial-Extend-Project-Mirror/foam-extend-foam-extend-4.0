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
    Foam::fixedDisplacementFvPatchVectorField

Description
    Time varying fixed displacement boundary condition with boundary
    non-orthogonal correction for stress solvers.
    Essentially fixedValue with non-orthogonal correction.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            timeVaryingFixedDisplacement;
        fileName        "$FOAM_CASE/time-series";
        outOfBounds     clamp;           // (error|warn|clamp|repeat)
    }
    @endverbatim

SourceFiles
    timeVaryingFixedDisplacementFvPatchVectorField.C

Authro
    Philip Cardiff
    Zeljko Tukovic

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingFixedDisplacementFvPatchVectorField_H
#define timeVaryingFixedDisplacementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedDisplacementFvPatchVectorField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class timeVaryingFixedDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingFixedDisplacementFvPatchVectorField
:
    public fixedDisplacementFvPatchVectorField
{

    // Private Data

        //- The time series being used, including the bounding treatment
        interpolationTable<vector> timeSeries_;

public:

    //- Runtime type information
    TypeName("timeVaryingFixedDisplacement");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingFixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingFixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // timeVaryingFixedDisplacementFvPatchVectorField
        //  onto a new patch
        timeVaryingFixedDisplacementFvPatchVectorField
        (
            const timeVaryingFixedDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingFixedDisplacementFvPatchVectorField
        (
            const timeVaryingFixedDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingFixedDisplacementFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingFixedDisplacementFvPatchVectorField
        (
            const timeVaryingFixedDisplacementFvPatchVectorField&,
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
                new timeVaryingFixedDisplacementFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return the time series used
            const interpolationTable<vector>& timeSeries() const
            {
                return timeSeries_;
            }

        // Evaluation functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
