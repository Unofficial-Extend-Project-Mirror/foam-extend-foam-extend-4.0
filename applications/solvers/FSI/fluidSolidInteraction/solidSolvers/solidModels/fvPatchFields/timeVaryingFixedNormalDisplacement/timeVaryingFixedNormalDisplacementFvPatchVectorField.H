/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    Foam::timeVaryingFixedNormalDisplacementFvPatchVectorField

Description
    Time varying fixed displacement boundary condition with boundary
    non-orthogonal correction for stress solvers.
    Essentially fixedValue with non-orthogonal correction.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            timeVaryingFixedNormalDisplacement;
        fileName        "$FOAM_CASE/time-series";
        outOfBounds     clamp;           // (error|warn|clamp|repeat)
    }
    @endverbatim

SourceFiles
    timeVaryingFixedNormalDisplacementFvPatchVectorField.C

Authro
    Philip Cardiff
    Zeljko Tukovic

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingFixedNormalDisplacementFvPatchVectorField_H
#define timeVaryingFixedNormalDisplacementFvPatchVectorField_H

#include "fixedNormalDisplacementFvPatchVectorField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class timeVaryingFixedNormalDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingFixedNormalDisplacementFvPatchVectorField
:
    public fixedNormalDisplacementFvPatchVectorField
{

    // Private Data

        //- The time series being used, including the bounding treatment
        interpolationTable<vector> timeSeries_;

public:

    //- Runtime type information
    TypeName("timeVaryingFixedNormalDisplacement");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingFixedNormalDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingFixedNormalDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // timeVaryingFixedNormalDisplacementFvPatchVectorField
        //  onto a new patch
        timeVaryingFixedNormalDisplacementFvPatchVectorField
        (
            const timeVaryingFixedNormalDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingFixedNormalDisplacementFvPatchVectorField
        (
            const timeVaryingFixedNormalDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingFixedNormalDisplacementFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingFixedNormalDisplacementFvPatchVectorField
        (
            const timeVaryingFixedNormalDisplacementFvPatchVectorField&,
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
                new timeVaryingFixedNormalDisplacementFvPatchVectorField
                (*this, iF)
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
