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
    Foam::fixedDisplacementFixedRotationFvPatchVectorField

Description
    time varying fixed displacement boundary condition and
    time varying fixed rotation boundary condition.

SourceFiles
    fixedDisplacementFixedRotationFvPatchVectorField.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef fixedDisplacementFixedRotationFvPatchVectorField_H
#define fixedDisplacementFixedRotationFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class fixedDisplacementFixedRotationFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedDisplacementFixedRotationFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

    // Private Data

        //- Name of the displacement field
        const word fieldName_;

        //- time varying displacement
        interpolationTable<vector> displacementTimeSeries_;

        //- Rotation in degrees
        //const scalar rotationAngle_;
        interpolationTable<scalar> rotationAngleTimeSeries_;

        //- Rotation axis
        const vector rotationAxis_;

        //- Origin of rotation
        const vector rotationOrigin_;


public:

    //- Runtime type information
    TypeName("fixedDisplacementFixedRotation");


    // Constructors

        //- Construct from patch and internal field
        fixedDisplacementFixedRotationFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedDisplacementFixedRotationFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // fixedDisplacementFixedRotationFvPatchVectorField
        // onto a new patch
        fixedDisplacementFixedRotationFvPatchVectorField
        (
            const fixedDisplacementFixedRotationFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedDisplacementFixedRotationFvPatchVectorField
        (
            const fixedDisplacementFixedRotationFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new fixedDisplacementFixedRotationFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedDisplacementFixedRotationFvPatchVectorField
        (
            const fixedDisplacementFixedRotationFvPatchVectorField&,
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
                new fixedDisplacementFixedRotationFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        // Access functions

            //- field name
            word fieldName() const
            {
                return fieldName_;
            }

        // Evaluation functions

            //- Return patch-normal gradient
            //- with non-orthogonal correction regardless of whether snGrad
            //- is corrected or not
            virtual tmp<Field<vector> > snGrad() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

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
