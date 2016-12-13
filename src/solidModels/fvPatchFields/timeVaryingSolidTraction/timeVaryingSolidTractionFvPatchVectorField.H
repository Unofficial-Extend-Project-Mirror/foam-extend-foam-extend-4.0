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
    timeVaryingSolidTractionFvPatchVectorField

Description
    Time-varying traction boundary condition for solid foam solvers.
    Suitable for for total strain or incremental strain total or
    updated Lagrangian approaches.
    Contains thermal-elastic, plastic, non-linear terms if required.
    And non-orthogonal correction.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            timeVaryingSolidTraction;
        field           U; // or DU
        nonLinear       OFF;
        fileName        "$FOAM_CASE/time-series";
        outOfBounds     clamp;           // (error|warn|clamp|repeat)
    }
    @endverbatim

Author
    Philip Cardiff

SourceFiles
    timeVaryingSolidTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingSolidTractionFvPatchVectorField_H
#define timeVaryingSolidTractionFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidTractionFvPatchVectorField.H"
#include "interpolationTable.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class timeVaryingSolidTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingSolidTractionFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{

    // Private Data

        //- The time series being used, including the bounding treatment
        interpolationTable<vector> timeSeries_;


public:

    //- Runtime type information
    TypeName("timeVaryingSolidTraction");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingSolidTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingSolidTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  timeVaryingSolidTractionFvPatchVectorField onto a new patch
        timeVaryingSolidTractionFvPatchVectorField
        (
            const timeVaryingSolidTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingSolidTractionFvPatchVectorField
        (
            const timeVaryingSolidTractionFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingSolidTractionFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingSolidTractionFvPatchVectorField
        (
            const timeVaryingSolidTractionFvPatchVectorField&,
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
                new timeVaryingSolidTractionFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the time series used
            const interpolationTable<vector>& timeSeries() const
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
