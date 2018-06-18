/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::blockSymmPlaneFvPatchVectorField

Description
    Block-coupled symmetry plane.  Experimental

Author
    Hrvoje Jasak, Wikki Ltd,  All rights reserved

SourceFiles
    blockSymmPlaneFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef blockSymmPlaneFvPatchVectorField_H
#define blockSymmPlaneFvPatchVectorField_H

#include "fvPatchFields.H"
#include "basicSymmetryFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class blockSymmPlaneFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class blockSymmPlaneFvPatchVectorField
:
    public basicSymmetryFvPatchVectorField
{
public:

    //- Runtime type information
    TypeName("blockSymmPlane");


    // Constructors

        //- Construct from patch and internal field
        blockSymmPlaneFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        blockSymmPlaneFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  blockSymmPlaneFvPatchVectorField onto a new patch
        blockSymmPlaneFvPatchVectorField
        (
            const blockSymmPlaneFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        blockSymmPlaneFvPatchVectorField
        (
            const blockSymmPlaneFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new blockSymmPlaneFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        blockSymmPlaneFvPatchVectorField
        (
            const blockSymmPlaneFvPatchVectorField&,
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
                new blockSymmPlaneFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field is block-coupled in the
            // boundary condition, ie. if the coupling coefficient is a
            // rank x rank implicit block
            virtual bool blockCoupled() const
            {
                return true;
            }


        // Evaluation functions

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            // Disabled.  HJ, 28/Mar/2016
            virtual tmp<vectorField> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            // Disabled.  HJ, 28/Mar/2016
            virtual tmp<vectorField> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            // Disabled.  HJ, 28/Mar/2016
            virtual tmp<vectorField> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            // Disabled.  HJ, 28/Mar/2016
            virtual tmp<vectorField> gradientBoundaryCoeffs() const;


        // Block-coupled coefficients

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<vectorCoeffField> blockValueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<vectorField> blockValueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<vectorCoeffField>
            blockGradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<vectorField> blockGradientBoundaryCoeffs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
