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
    Foam::fixedRotationFvPatchVectorField

Description
    Fixed rotation boundary condition with boundary non-orthogonal
    correction for stress solvers.

SourceFiles
    fixedRotationFvPatchVectorField.C

Author
    Philip Cardiff

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingFixedRotationFvPatchVectorField_H
#define timeVaryingFixedRotationFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "interpolationTable.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class timeVaryingFixedRotationFvPatch Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingFixedRotationFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

    // Private Data

        interpolationTable<scalar> timeSeries_;

        //- Name of the displacement field
        const word fieldName_;

        //- total accumulated angle - only used for UL
        scalar accumulatedAngle_;

        //- initial patch face centres
        const vectorField origFaceCentres_;

        //- Rotation in degrees
        //const scalar rotationAngle_;

        //- Rotation in degrees
        const vector rotationAxis_;

        //- Origin of rotation
        const vector rotationOrigin_;



public:

    //- Runtime type information
    TypeName("timeVaryingFixedRotation");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingFixedRotationFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingFixedRotationFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // timeVaryingFixedRotationFvPatchVectorField
        //  onto a new patch
        timeVaryingFixedRotationFvPatchVectorField
        (
            const timeVaryingFixedRotationFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingFixedRotationFvPatchVectorField
        (
            const timeVaryingFixedRotationFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingFixedRotationFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingFixedRotationFvPatchVectorField
        (
            const timeVaryingFixedRotationFvPatchVectorField&,
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
                new timeVaryingFixedRotationFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

              //- Return the time series used
              const interpolationTable<scalar>& timeSeries() const
              {
                  return timeSeries_;
              }


        // Evaluation functions

            //- Return patch-normal gradient
            //- with non-orthogonal correction regardless of whether snGrad
            //- is corrected or not
            virtual tmp<Field<vector> > snGrad() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
