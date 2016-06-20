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

#ifndef fixedRotationFvPatchVectorField_H
#define fixedRotationFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "nonLinearGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class fixedRotationFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedRotationFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

    // Private Data

        //- Name of the displacement field
        const word fieldName_;

        //- Rotation in degrees
        const scalar rotationAngle_;

        //- Rotation axis
        const vector rotationAxis_;

        //- Origin of rotation
        const vector rotationOrigin_;

        //- Initial mesh position
        const vectorField origCf_;

        //- Is it a non linear solver
        nonLinearGeometry::nonLinearType nonLinear_;


public:

    //- Runtime type information
    TypeName("fixedRotation");


    // Constructors

        //- Construct from patch and internal field
        fixedRotationFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedRotationFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedRotationFvPatchVectorField
        //  onto a new patch
        fixedRotationFvPatchVectorField
        (
            const fixedRotationFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedRotationFvPatchVectorField
        (
            const fixedRotationFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new fixedRotationFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedRotationFvPatchVectorField
        (
            const fixedRotationFvPatchVectorField&,
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
                new fixedRotationFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        // Access functions

            const nonLinearGeometry::nonLinearType& nonLinear() const
            {
                return nonLinear_;
            }

            nonLinearGeometry::nonLinearType& nonLinear()
            {
                return nonLinear_;
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
