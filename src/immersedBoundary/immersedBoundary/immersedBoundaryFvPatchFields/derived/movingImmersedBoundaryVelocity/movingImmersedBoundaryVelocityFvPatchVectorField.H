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
    Foam::movingImmersedBoundaryVelocityFvPatchVectorField

Description
    Foam::movingImmersedBoundaryVelocityFvPatchVectorField

SourceFiles
    movingImmersedBoundaryVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef movingImmersedBoundaryVelocityFvPatchVectorField_H
#define movingImmersedBoundaryVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "immersedBoundaryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class movingImmersedBoundaryVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class movingImmersedBoundaryVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Local reference to immersed boundary patch
        const immersedBoundaryFvPatch& ibPatch_;


        // Dead cell controls

            //- Set dead cell value
            Switch setDeadValue_;

            //- Dead cell value
            vector deadValue_;


    // Private Member Functions

        //- Set value in dead cells
        void setDeadValues();


public:

    //- Runtime type information
    TypeName("movingImmersedBoundaryVelocity");


    // Constructors

        //- Construct from patch and internal field
        movingImmersedBoundaryVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        movingImmersedBoundaryVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  movingImmersedBoundaryVelocityFvPatchVectorField
        //  onto a new patch
        movingImmersedBoundaryVelocityFvPatchVectorField
        (
            const movingImmersedBoundaryVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        movingImmersedBoundaryVelocityFvPatchVectorField
        (
            const movingImmersedBoundaryVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new movingImmersedBoundaryVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        movingImmersedBoundaryVelocityFvPatchVectorField
        (
            const movingImmersedBoundaryVelocityFvPatchVectorField&,
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
                new movingImmersedBoundaryVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return reference to immersed boundary patch
            const immersedBoundaryFvPatch& ibPatch() const
            {
                return ibPatch_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
