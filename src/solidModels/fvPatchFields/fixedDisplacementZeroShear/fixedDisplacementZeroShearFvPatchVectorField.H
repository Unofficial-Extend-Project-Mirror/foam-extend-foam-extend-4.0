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
    fixedDisplacementZeroShearFvPatchVectorField

Description
    The component of the displacement normal to the patch is applied,
    and the shear traction is set to zero.
    Patch contains non-orthogonal correction.

SourceFiles
    fixedDisplacementZeroShearFvPatchVectorField.C

Author
    Philip Cardiff
    Based on velocityDisplacement by Aleksandar Karac

\*---------------------------------------------------------------------------*/

#ifndef fixedDisplacementZeroShearFvPatchVectorField_H
#define fixedDisplacementZeroShearFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "nonLinearGeometry.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedDisplacementZeroShearFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{
    // Private data

        //- Name of the displacement field
        const word fieldName_;

        //- Is it a non linear solver
        nonLinearGeometry::nonLinearType nonLinear_;

        //- if it is an orthropic solver
        Switch orthotropic_;

public:

    //- Runtime type information
    TypeName("fixedDisplacementZeroShear");


    // Constructors

        //- Construct from patch and internal field
        fixedDisplacementZeroShearFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedDisplacementZeroShearFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // fixedDisplacementZeroShearFvPatchVectorField
        //  onto a new patch
        fixedDisplacementZeroShearFvPatchVectorField
        (
            const fixedDisplacementZeroShearFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new fixedDisplacementZeroShearFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedDisplacementZeroShearFvPatchVectorField
        (
            const fixedDisplacementZeroShearFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new fixedDisplacementZeroShearFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~fixedDisplacementZeroShearFvPatchVectorField()
        {}


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
