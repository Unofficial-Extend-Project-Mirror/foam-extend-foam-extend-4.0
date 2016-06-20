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
    Foam::elasticWallVelocityFvPatchVectorField

Description
    Foam::elasticWallVelocityFvPatchVectorField

SourceFiles
    elasticWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef elasticWallVelocityFvPatchVectorField_H
#define elasticWallVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class elasticWallVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class elasticWallVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

    //To keep track when to copy patch face centres
    label myTimeIndex_;

    //Old face centres of patch. Needed for (second order) time integration
    vectorField Fc_;
    vectorField oldFc_;
    vectorField oldOldFc_;


public:

    //- Runtime type information
    TypeName("elasticWallVelocity");


    // Constructors

        //- Construct from patch and internal field
        elasticWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        elasticWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given elasticWallVelocityFvPatchVectorField
        //  onto a new patch
        elasticWallVelocityFvPatchVectorField
        (
            const elasticWallVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        elasticWallVelocityFvPatchVectorField
        (
            const elasticWallVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new elasticWallVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        elasticWallVelocityFvPatchVectorField
        (
            const elasticWallVelocityFvPatchVectorField&,
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
                new elasticWallVelocityFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

    // Evaluation functions

        //- Return patch-normal gradient
        virtual tmp<Field<vector> > snGrad() const;

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
